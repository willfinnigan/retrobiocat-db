

# delete all activity data except for finnigan paper
# delete all sequences
# delete all papers
# delete all enzyme types except CAR
# delete all reaction rules except carboxylic acid reduction
# delete all users
# delete all issues + suggestions ect

from retrobiocat_web.mongo.model_queries.activity_queries import activity_in_paper
from retrobiocat_web.mongo.model_queries.paper_queries import paper_from_doi
from retrobiocat_web.mongo.models.biocatdb_models import Activity, Sequence, Paper, Molecule, MoleculeStats, Tag, \
    UniRef50, SSN_record, ActivityIssue, ActivityMol, IrrelevantPaper, EnzymeTypeSearchHistory
from retrobiocat_web.mongo.models.comments import Comment
from retrobiocat_web.mongo.models.enzyme_type import EnzymeType
from retrobiocat_web.mongo.models.other_models import AutoJobStatus, Chebi, Rhea, PaperSuggestion
from retrobiocat_web.mongo.models.reaction_misc_models import Issue
from retrobiocat_web.mongo.models.reaction_models import Reaction
from retrobiocat_web.mongo.models.retro_tests import TestCascade, TestCascadeResult, TestCascadeRun
from retrobiocat_web.mongo.models.user_models import User, Role
from retrobiocat_web.mongo.models.user_saves import Network, MyMolecule


def delete_all_data_except_paper(paper_to_keep):
    activities = Activity.objects(paper__ne=paper_to_keep)
    for act in activities:
        act.delete()

    sequences = Sequence.objects(papers__ne=paper_to_keep)
    for seq in sequences:
        seq.delete()

    papers = Paper.objects(doi__ne=paper_to_keep.doi)
    for paper in papers:
        paper.delete()

    print('Done')

def delete_all_enzyme_types_except(et_abbrev_to_keep):
    enzyme_types = EnzymeType.objects(enzyme_type__ne=et_abbrev_to_keep)
    for et in enzyme_types:
        et.delete()

    print('Done')

def delete_all_rules_except(reaction_name_to_keep):
    reactions = Reaction.objects(name__ne=reaction_name_to_keep)
    for r in reactions:
        r.delete()

    print('Done')

def delete_all_misc_data():
    MoleculeStats.drop_collection()
    Tag.drop_collection()
    UniRef50.drop_collection()
    SSN_record.drop_collection()
    ActivityIssue.drop_collection()
    ActivityMol.drop_collection()
    IrrelevantPaper.drop_collection()
    EnzymeTypeSearchHistory.drop_collection()
    Comment.drop_collection()
    AutoJobStatus.drop_collection()
    Chebi.drop_collection()
    Rhea.drop_collection()
    PaperSuggestion.drop_collection()
    Issue.drop_collection()
    TestCascade.drop_collection()
    TestCascadeResult.drop_collection()
    TestCascadeRun.drop_collection()
    Network.drop_collection()
    MyMolecule.drop_collection()

    print('done')

def drop_all_users():
    User.drop_collection()
    Role.drop_collection()

def drop_test_cascases():
    TestCascade.drop_collection()
    TestCascadeResult.drop_collection()
    TestCascadeRun.drop_collection()


def set_owners_to_none():
    for paper in Paper.objects():
        paper.owner = None
        paper.added_by = None
        paper.edits_by = []
        paper.seq_reviewed_by = None
        paper.activity_reviewed_by = None
        paper.save()

    for seq in Sequence.objects():
        seq.owner = None
        seq.added_by = None
        seq.edits_by = []
        seq.save()

    for act in Activity.objects():
        act.added_by = None
        act.edits_by = None






if __name__ == '__main__':
    from retrobiocat_web.mongo.default_connection import make_default_connection
    make_default_connection()

    doi = '10.1002/cctc.201601249'
    paper = paper_from_doi(doi)

    delete_all_data_except_paper(paper)
    delete_all_enzyme_types_except('CAR')
    delete_all_rules_except('Carboxylic acid reduction')
    delete_all_misc_data()
    drop_all_users()
    drop_test_cascases()
    set_owners_to_none()
