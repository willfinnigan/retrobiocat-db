
class Visualiser_Config():

    def __init__(self):
        self.molSize = (300,300)
        self.display_cofactors = True
        self.colour_arrows = "None"  # or Complexity change
        self.colour_reactions = "Substrate specificity"  # no other options currently
        self.show_negative_enzymes = True
        self.fix_target = False

    def update_from_dict(self, attr_dict):
        current_dict = self.to_dict()
        for key, value in attr_dict.items():
            if key in current_dict:
                setattr(self, key, value)
        return self

    def to_dict(self):
        return self.__dict__

