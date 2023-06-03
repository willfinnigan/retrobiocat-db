import mongoengine as db
from datetime import datetime


class Reaction(db.Document):
    name = db.StringField(unique=True)
    smarts = db.ListField(db.StringField())
    enzyme_types = db.ListField(db.StringField())
    cofactors = db.DictField()
    positive_tests = db.ListField(db.StringField())
    negative_tests = db.ListField(db.StringField())
    example_rxn_string = db.StringField()
    type = db.StringField()
    experimental = db.BooleanField(default=False)
    two_step = db.BooleanField(default=False)
    requires_absence_of_water = db.BooleanField(default=False)
    steps = db.ListField(db.ListField(db.ReferenceField('self')))
    last_updated = db.DateTimeField(default=datetime.utcnow)
    product_seeds = db.ListField(db.StringField())
    substrate_1_seeds = db.ListField(db.StringField())
    substrate_2_seeds = db.ListField(db.StringField())


if __name__ == '__main__':
    from retrobiocat_web.mongo.default_connection import make_default_connection
    make_default_connection(host='138.68.135.53', database='test')

    #redam = Reaction.objects(name='Reductive amination')[0]
    #print(redam.to_json())

    imine_formation = Reaction.objects(name='Imine formation')[0]
    imine_reduction = Reaction.objects(name='Imine reduction')[0]

    iminium_formation = Reaction.objects(name='Iminium formation')[0]
    iminium_reduction = Reaction.objects(name='Iminium reduction')[0]

    redam = Reaction.objects(name='Reductive amination').select_related()[0]

    print(redam.to_json())

    redam.two_step = True
    redam.experimental = False
    redam.smarts = []
    redam.steps = [[imine_reduction, imine_formation], [iminium_reduction, iminium_formation]]
    redam.save()




