import numpy as np
from bokeh.models import ColumnDataSource, HoverTool, CustomJS, TapTool
from bokeh.plotting import figure
from retrobiocat_web.analysis.data_query import get_data
from retrobiocat_web.logging import add_logger


class Temperature_Processor():

    def __init__(self, room_temp=25, log_level='WARNING'):
        self.room_temp = room_temp
        self.logger = add_logger('TemperatureProcessor', level=log_level)

    def parse_temps(self, temperatures):
        """
        Parse the temperatures so that all numbers
        - remove empty entries ('')
        - RT should be 25
        - all floats
        """

        parsed_temps = []
        for t in temperatures:

            t = self.parse_single_temp(t)

            if t is not None:
                parsed_temps.append(t)

        return parsed_temps

    def parse_single_temp(self, temp):
        t = self._temp_to_float(temp)

        if isinstance(t, str):
            t = self._remove_spaces(t)
            t = self._remove_degrees_and_convert_to_int(t)
            t = self._room_temp_to_number(t)

        t = self._nan_to_none(t)

        if isinstance(t, float):
            return t
        else:
            return None

    def _nan_to_none(self, t):
        if np.isnan(t):
            return None
        return t



    def _room_temp_to_number(self, t):
        """Converts a string representing room temperature to the default RT number"""

        # if any of these strings are present do conversion
        room_temp_strings = ['room', 'rt', 'r.t', 'r.t.', 'ambient']

        if isinstance(t, str):
            if t.lower() in room_temp_strings:
                t = self.room_temp
            elif t == '':
                pass
            else:
                self.logger.debug(f'{t} is a string but doesnt match rt strings')
        return t

    def _remove_spaces(self, t):
        return t.replace(' ', '')

    def _remove_degrees_and_convert_to_int(self, t):
        def _remove_deg_text(t, text):
            if text in t:
                t = t.replace(text, '')
            return t

        if isinstance(t, str):
            t = _remove_deg_text(t, '℃')
            t = _remove_deg_text(t, '°C')
            t = self._temp_to_float(t)

        return t

    def _temp_to_float(self, t):
        """Converts t to an int or returns a None"""
        try:
            return float(t)
        except:
            return t


class TemperatureHistrogramCreator():

    def __init__(self):
        self.temperature_processor = Temperature_Processor()

    def _get_temps(self, df):
        if 'temperature' not in df.columns:
            return []
        return list(df.temperature)

    def create_graph(self, data_query, bins=15):
        # histrogram of temperatures
        # clicking on a bar will launch...?  activity data in that group.

        df = data_query.get_activity_df()
        temperatures = self._get_temps(df)
        try:
            parsed_temperatures = [t for t in temperatures if not np.isnan(t)]
        except:
            print('ERROR PARSING TEMPERATURES')
            parsed_temperatures = []

        hist, edges = np.histogram(parsed_temperatures, bins=bins)

        data = {'counts': hist,
                'lefts': edges[:-1],
                'rights': edges[1:]}
        source = ColumnDataSource(data)

        p = figure(title=f'Reaction temperatures histogram ({len(parsed_temperatures)} out of {len(temperatures)} rows with valid data)', tools='', toolbar_location=None,
                   frame_width=500, frame_height=400)
        p.quad(source=source, top='counts', bottom=0, left='lefts', right='rights',
               fill_color="#b00000", line_width=1.5, line_color='black', name='temp_hist',
               hover_alpha=0.65)

        p.xaxis.axis_label = r"$$Temperature \degree C$$"
        p.yaxis.axis_label = r'$$Count$$'
        p.y_range.start = 0
        p.xgrid.visible = False
        p.ygrid.visible = False
        p.outline_line_width = 1
        p.outline_line_color = "black"

        HT = HoverTool(names=['temp_hist'])
        HT.tooltips = [('Count', '@counts'),
                       ('Temperature', '@lefts to @rights')]
        HT.mode = 'mouse'
        p.add_tools(HT)

        callback = CustomJS(args=dict(source=source),
                            code=f"""var ind = cb_data.source.selected.indices[0];
                                             var temp_left = Math.round(source.data.lefts[ind] * 100) / 100
                                             var temp_right = Math.round(source.data.rights[ind] * 100) / 100
                                             launch_temp_modal(temp_left, temp_right)""")

        TT = TapTool(names=['temp_hist'])
        TT.callback = callback
        p.add_tools(TT)
        return p


if __name__ == '__main__':
    from retrobiocat_web.mongo.default_connection import make_default_connection

    make_default_connection()
    from bokeh.plotting import output_file, save

    dq = get_data.DataQuery(enzyme_type='CAR', reaction='Carboxylic acid reduction', log_level=1)
    temp_hist_maker = TemperatureHistrogramCreator()
    temp_graph = temp_hist_maker.create_graph(dq)

    output_file('temp_graph_test.html')
    save(temp_graph)




