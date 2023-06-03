import requests
import pandas as pd
import io
from retrobiocat_web.mongo.models.biocatdb_models import UniRef50, EnzymeType
import copy
import time

def uni100_for_enzyme_type(enzyme_type_obj):
    unirefs = UniRef50.objects(enzyme_type=enzyme_type_obj)
    code_list = []
    for ur in unirefs:
        code_list += ur.uni100
    return code_list

def uni100_for_specified_uniref50(list_ur):
    unirefs = UniRef50.objects(enzyme_name__in=list_ur)
    code_list = []
    for ur in unirefs:
        code_list += ur.uni100
    return code_list

class BulkUniProt(object):

    def __init__(self, log_level=0):
        self.url = 'https://www.uniprot.org/uploadlists/'
        self.params = {'from': 'ACC+ID',
                       'to': 'ACC',
                       'format': 'tab',
                       'columns': 'sequence',
                       'query': ''}

        self.batch_size = 500
        self.log_level = log_level

    # Public methods
    def get_seqs(self, code_list):
        chunked_code_list = self._divide_chunks(code_list)

        list_dfs = []
        for cl in chunked_code_list:
            text = self._uniprot_request(cl)
            list_dfs.append(self._response_text_to_df(text))

        df = pd.concat(list_dfs)

        return df

    # Private methods
    def _uniprot_request(self, query_list, retry_num=3, retry_wait=10):
        params = copy.deepcopy(self.params)
        params['query'] = " ".join(query_list)

        for i in range(retry_num):
            try:
                response = requests.get(self.url, params=params)

                if response.status_code in [200]:
                    return response.text
                else:
                    self._log(f"UniProt bulk download error - Status code = {response.status_code}", level=0)
                    return ''

            except:
                time.sleep(retry_wait)

        self._log(f"UniProt bulk download error - time out after {retry_num} retries with {retry_wait} seconds between", level=0)
        return ''

    @staticmethod
    def _response_text_to_df(text):
        if text == "":
            return pd.DataFrame()

        return pd.read_csv(io.StringIO(text), sep='\t', names=['sequence', 'id'], header=0)

    def _divide_chunks(self, list_to_chunk):
        for i in range(0, len(list_to_chunk), self.batch_size):
            yield list_to_chunk[i:i + self.batch_size]

    def _log(self, msg, level=1):
        if level <= self.log_level:
            print(f"UniProt_Bulk_Download: {msg}")

if __name__ == '__main__':
    from retrobiocat_web.mongo.default_connection import make_default_connection
    make_default_connection(host='138.68.135.53', database='test')
    enzyme_type = EnzymeType.objects(enzyme_type='CAR')[0]

    # this generates a list of 10 id's to test work flow on
    unirefs = UniRef50.objects(enzyme_type=enzyme_type).distinct('enzyme_name')[0:10]

    code_list = uni100_for_specified_uniref50(unirefs)

    df = BulkUniProt().get_seqs(code_list)

    print(df.head())
    print(len(df.index))

    #df = pd.read_csv(io.StringIO(response.text), sep='\t', names=['sequence',  'id'], header=0)
    #print(df.head())

