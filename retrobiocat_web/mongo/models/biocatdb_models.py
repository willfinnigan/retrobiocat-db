import mongoengine as db
from rdkit import Chem
from datetime import datetime

from retrobiocat_web.mongo.models.enzyme_type import EnzymeType
from retrobiocat_web.mongo.models.user_models import User
from retrobiocat_web.mongo.models.comments import Comment


class Paper(db.Document):
    doi = db.StringField(unique=True)
    short_citation = db.StringField(default='')
    html = db.StringField(default='')
    owner = db.ReferenceField(User, reverse_delete_rule=1)
    added_by = db.ReferenceField(User, reverse_delete_rule=1)
    edits_by = db.ListField(db.ReferenceField(User, reverse_delete_rule=4))

    title = db.StringField(default='')
    authors = db.ListField(db.StringField())
    journal = db.StringField(default='')
    date = db.DateField()

    status = db.StringField(default='')
    tags = db.ListField(db.StringField())

    metadata_reviewed = db.BooleanField(default=False)

    seq_review_ready = db.BooleanField(default=False)
    seq_reviewed = db.BooleanField(default=False)
    seq_reviewed_by = db.ReferenceField(User, reverse_delete_rule=1)

    activity_review_ready = db.BooleanField(default=False)
    activity_reviewed = db.BooleanField(default=False)
    activity_reviewed_by = db.ReferenceField(User, reverse_delete_rule=1)

    # if the paper has been added by a standard user, is it reviewed for inclusion in database?
    reviewed = db.BooleanField(default=True)

    # legacy field, no longer needed
    reviewed_by = db.ReferenceField(User, reverse_delete_rule=1)

    has_issues = db.BooleanField(default=False)
    comments = db.ListField(db.ReferenceField(Comment, reverse_delete_rule=4))
    high_importance = db.BooleanField(default=False)
    no_relevant_activity = db.BooleanField(default=False)

    def authors_str(self):
        return ", ".join(str(x) for x in self.authors)

    def owner_text(self):
        self.select_related()
        if self.owner is None:
            return ''
        return f"{self.owner.first_name} {self.owner.last_name}, {self.owner.affiliation}"

    def tags_str(self):
        return ", ".join(str(x) for x in self.tags)

    def v_short_cit(self):
        if len(self.authors) != 0:
            return f"{self.authors[0]} et al"
        return self.title[0:10] + '..'

    def str_date(self):
        if self.date != None:
            return self.date.strftime('%B-%Y')
        return None

    def __unicode__(self):
        return self.short_citation

    def __str__(self):
        return self.short_citation

class SequenceOtherNamesData(db.EmbeddedDocument):
    name = db.StringField()
    n_tag = db.StringField()
    c_tag = db.StringField()
    notes = db.StringField()
    papers = db.ListField(db.ReferenceField(Paper))

class Sequence(db.Document):
    enzyme_type = db.StringField(max_length=120, required=True)
    enzyme_name = db.StringField(max_length=120, unique=True, required=True)
    other_names = db.ListField(db.StringField())  # not used any more, use other_names_data instead
    other_names_data = db.EmbeddedDocumentListField(SequenceOtherNamesData, default=[])
    n_tag = db.StringField()
    sequence = db.StringField(default='')
    seqvec = db.StringField()
    unirep = db.StringField()
    c_tag = db.StringField()
    sequence_unavailable = db.BooleanField(default=False)

    accession = db.StringField(max_length=20, default='')
    other_identifiers = db.ListField(db.StringField(max_length=20))

    structure = db.BooleanField(default=False)
    pdb = db.StringField(default='')
    mutant_of = db.StringField(default='')
    notes = db.StringField(default='')
    bioinformatics_ignore = db.BooleanField(default=False)

    owner = db.ReferenceField(User, reverse_delete_rule=1)
    added_by = db.ReferenceField(User, reverse_delete_rule=1)
    edits_by = db.ListField(db.ReferenceField(User, reverse_delete_rule=4))
    papers = db.ListField(db.ReferenceField(Paper, reverse_delete_rule=4))

    blast = db.DateTimeField(default=None)
    alignments_made = db.BooleanField(default=False)

    objects_to_update = []
    reviewed = db.BooleanField(default=False)

    def __unicode__(self):
        return self.enzyme_name

    def __str__(self):
        return self.enzyme_name

    def owner_text(self):
        self.select_related()
        if self.owner is None:
            return ''
        return f"{self.owner.first_name} {self.owner.last_name}, {self.owner.affiliation}"

    meta = {'indexes': ['enzyme_name', 'enzyme_type']}


class Molecule(db.DynamicDocument):
    smiles = db.StringField()
    mol = db.BinaryField()
    #mfp2_2048 = db.StringField()
    #rdfp2_2048 = db.StringField()

    def get_mol(self):
        return Chem.Mol(self.mol)

    meta = {'indexes': ['smiles']}

class MoleculeStats(db.Document):
    smiles = db.StringField()
    descriptors = db.DictField()


class Activity(db.Document):
    enzyme_type = db.StringField()
    enzyme_name = db.StringField()
    reaction = db.StringField()
    short_citation = db.StringField()
    html_doi = db.StringField()

    added_by = db.ReferenceField(User, reverse_delete_rule=1)
    added_by_string = db.StringField()  # legacy field, no longer in use
    paper = db.ReferenceField(Paper, reverse_delete_rule=1)
    edits_by = db.ListField(db.ReferenceField(User, reverse_delete_rule=4))

    substrate_1_smiles = db.StringField()
    substrate_2_smiles = db.StringField()
    product_1_smiles = db.StringField()
    cascade_num = db.StringField()

    temperature = db.FloatField()
    ph = db.FloatField()
    solvent = db.StringField()
    other_conditions = db.StringField()
    notes = db.StringField()
    reaction_vol = db.StringField()
    formulation = db.StringField()
    biocat_conc = db.StringField()

    kcat = db.FloatField()
    km = db.FloatField()
    mw = db.FloatField()

    substrate_1_conc = db.StringField()
    substrate_2_conc = db.StringField()
    specific_activity = db.FloatField()
    conversion = db.FloatField()
    conversion_time = db.FloatField()

    categorical = db.StringField()
    binary = db.BooleanField()

    selectivity = db.StringField()
    auto_generated = db.BooleanField()
    reviewed = db.BooleanField(default=False)

    meta = {'indexes': ['enzyme_name', 'enzyme_type', 'reaction']}

class Tag(db.Document):
    seq = db.StringField()
    n_term = db.BooleanField(default=False)
    c_term = db.BooleanField(default=False)

    def __unicode__(self):
        return self.seq

    def __str__(self):
        return self.seq

class UniRef50(db.Document):
    # These fields are shared with Sequence
    enzyme_name = db.StringField()
    sequence = db.StringField()

    enzyme_type = db.ReferenceField(EnzymeType, reverse_delete_rule=2)
    result_of_blasts_for = db.ListField(db.ReferenceField(Sequence, reverse_delete_rule=4))
    seq_match = db.ReferenceField(Sequence, reverse_delete_rule=1)
    blast_round = db.IntField()
    alignments_made = db.BooleanField(default=False)

    protein_name = db.StringField()
    tax = db.StringField()
    tax_id = db.StringField()

    num_uni90 = db.IntField()
    num_uni100 = db.IntField()
    num_uniprot = db.IntField()
    uni100 = db.ListField(db.StringField())
    cluster_name = db.StringField()
    rep_name = db.StringField()
    pfam_codes = db.ListField(db.StringField())
    pfams = db.DictField()
    sp_annotated = db.BooleanField()
    pdbs = db.ListField(db.StringField())
    rhea = db.ListField(db.StringField())

    meta = {'indexes': ['enzyme_type', 'enzyme_name']}

    def __unicode__(self):
        return self.enzyme_name

    def __str__(self):
        return self.enzyme_name

class UniRef90(db.DynamicDocument):
    pass

class Alignment(db.DynamicDocument):
    pass

class SeqSimNet(db.DynamicDocument):
    pass

class SSN_record(db.Document):
    enzyme_type = db.ReferenceField(EnzymeType, reverse_delete_rule=2)
    status = db.StringField(default='Newly created')

    precalc_status = db.StringField(default='Empty')
    precalculated_vis = db.DictField(default={})
    num_at_alignment_score = db.DictField(default={})
    identity_at_alignment_score = db.DictField(default={})
    pos_at_alignment_score = db.DictField(default={})

    meta = {'indexes': ['enzyme_type']}

class ActivityIssue(db.Document):
    activity = db.ReferenceField(Activity)
    raised_by = db.ReferenceField(User, reverse_delete_rule=2)
    status = db.StringField(default='Open')
    comments = db.ListField(db.ReferenceField(Comment, reverse_delete_rule=4))
    date = db.DateTimeField(default=datetime.utcnow)

class ActivityMol(db.Document):
    paper = db.ReferenceField(Paper)
    name = db.StringField(default="")
    chem_name = db.StringField(default="")
    smi = db.StringField()
    svg = db.StringField()

class IrrelevantPaper(db.Document):
    doi = db.StringField(required=True)
    enzyme_type = db.ReferenceField(EnzymeType, reverse_delete_rule=2)

class EnzymeTypeSearchHistory(db.Document):
    enzyme_type = db.ReferenceField(EnzymeType, reverse_delete_rule=2)
    date = db.DateTimeField(default=datetime.utcnow())
    user = db.ReferenceField(User, reverse_delete_rule=1)
    score = db.IntField()


if __name__ == '__main__':
    from retrobiocat_web.mongo.default_connection import make_default_connection
    make_default_connection(host='138.68.135.53', database='test')



