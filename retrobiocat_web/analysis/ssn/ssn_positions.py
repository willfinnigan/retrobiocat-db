
class ClusterPositioner(object):

    def __init__(self):

        self.v_move = 0
        self.max_width = 10000
        self.gutter = 500
        self.current_location = [0, 0]

    def get_cluster_dimensions(self, pos_dict):
        x = []
        y = []
        for node, positions in pos_dict.items():
            x.append(positions[0])
            y.append(positions[1])

        max_x = max(x)
        max_y = max(y)
        min_x = min(x)
        min_y = min(y)

        x_dim = max_x - min_x
        y_dim = max_y - min_y
        center = (max_x-(x_dim/2), max_y-(y_dim/2))

        return x_dim, y_dim, center

    @staticmethod
    def move(pos_dict, move):
        new_pos_dict = {}
        for key in pos_dict:
            new_pos_dict[key] = [pos_dict[key][0] + move[0], pos_dict[key][1] + move[1]]
        return new_pos_dict

    @staticmethod
    def make_coords_positive(pos_dict):
        new_pos_dict = {}
        move_x, move_y = 0, 0
        for key, value in pos_dict.items():
            x, y = value[0], value[1]
            if x > move_x:
                move_x = x
            if y > move_y:
                move_y = y

        for key, value in pos_dict.items():
            x, y = value[0], value[1]
            x += move_x
            y += move_y
            new_pos_dict[key] = [x, y]

        return new_pos_dict

    def position(self, pos_dict):
        pos_dict = self.make_coords_positive(pos_dict)

        h_dim, v_dim, center = self.get_cluster_dimensions(pos_dict)

        # should this be a new row?
        if self.current_location[0] + h_dim > self.max_width:
            self.current_location[0] = 0
            self.current_location[1] += self.v_move + self.gutter + v_dim/2
            self.v_move = 0

        # center nodes around zero
        pos_dict = self.move(pos_dict, (-center[0], -center[1]))

        self.current_location[0] += h_dim/2
        pos_dict = self.move(pos_dict, self.current_location)

        # add space between clusters
        self.current_location[0] += h_dim/2
        self.current_location[0] += self.gutter

        if v_dim/2 > self.v_move:
            self.v_move = v_dim/2

        return pos_dict

    def round_positions(self, pos_dict, round_to=2):
        rounded_pos_dict = {}
        for key in pos_dict:
            rounded_pos_dict[key] = [round(pos_dict[key][0], round_to),
                                     round(pos_dict[key][1], round_to)]
        return rounded_pos_dict