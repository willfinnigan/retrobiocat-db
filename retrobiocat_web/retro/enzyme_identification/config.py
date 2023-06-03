

class Specificity_Scorer_Config():

    def __init__(self):
        self.run_similarity_search = True
        self.similarity_threshold = 0.55
        self.products_only = True
        self.only_active = True
        self.only_reviewed = True

    def update_from_dict(self, attr_dict):
        current_dict = self.to_dict()
        for key, value in attr_dict.items():
            if key in current_dict:
                setattr(self, key, value)
        return self

    def to_dict(self):
        return self.__dict__

