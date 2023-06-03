from retrobiocat_web.mongo.modal_updates import bioinformatics_CUD
from retrobiocat_web.mongo.model_queries import reaction_queries, sequence_queries, activity_queries
from retrobiocat_web.mongo.models.biocatdb_models import Sequence, Activity, Paper
from retrobiocat_web.analysis.ssn import ssn_main
from retrobiocat_web.mongo.models.enzyme_type import EnzymeType
from retrobiocat_web.mongo.models.reaction_models import Reaction


def relabel_enzyme_type_data(enz_type, new_name, rename_ssn=True):
    """Rename all enzyme_type embedded fields"""
    try:
        old_name = enz_type.enzyme_type
        if rename_ssn == True:
            bioinformatics_CUD.rename_ssn(old_name, new_name)

        sequences = Sequence.objects(enzyme_type=old_name)
        for seq in sequences:
            seq.enzyme_type = new_name
            seq.save()

        # if there were sequences renamed, then bioinformatics will need updating
        if len(sequences) != 0:
            bioinformatics_CUD.enzyme_type_needs_blasting(enz_type.enzyme_type)
            bioinformatics_CUD.SSN_needs_updating(enz_type.enzyme_type)

        for reaction in Reaction.objects(enzyme_types=old_name):
            reaction.enzyme_types.remove(old_name)
            reaction.enzyme_types.append(new_name)
            reaction.cofactors[new_name] = reaction.cofactors.pop(old_name)
            reaction.save()
        for activity in Activity.objects(enzyme_type=old_name):
            activity.enzyme_type = new_name
            activity.save()
        for paper in Paper.objects(tags=old_name):
            for i, tag in enumerate(paper.tags):
                if tag == old_name:
                    paper.tags[i] = new_name
            paper.save()
        return []
    except Exception as e:
        return ['Error renaming enzyme_type embedded fields', str(e)]

def change_enzyme_type_name(enz_type, new_name, rename_ssn=True):
    """
    Renaming an enzyme type means this field must be updated in Sequence, Reaction and Activity also
    Paper tags must also be renamed
    """
    old_name = enz_type.enzyme_type

    if new_name in list(EnzymeType.objects().distinct('enzyme_type')):
        return [f'Enzyme type with the abbreviation {new_name} already exists']

    # relabel all the embedded data
    issues = relabel_enzyme_type_data(enz_type, new_name, rename_ssn=rename_ssn)
    if len(issues) != 0:
        return issues

    # if relabelling was ok, then rename the enzyme_type
    try:
        enz_type.enzyme_type = new_name
        enz_type.save()
        return []
    except Exception as e:
        return [f'Error renaming enzyme type from {old_name} to {new_name}', str(e)]


def update_enzyme_type(enz_type, new_name, description, full_name, other_abbreviations, rep_reaction_name):
    """Update the fields for an enzyme_type"""

    try:
        issues = []
        if other_abbreviations == '':
            other_abbreviations = None
        else:
            other_abbreviations = other_abbreviations.split(', ')

        # get reaction object for representative reaction
        rep_reaction = reaction_queries.reaction_from_name(rep_reaction_name)  # will be none if not found

        # update fields
        enz_type.description = description
        enz_type.rep_reaction = rep_reaction
        enz_type.full_name = full_name
        enz_type.other_abbreviations = other_abbreviations
        enz_type.save()

        # update the enzyme_type abbreviation if its not identical
        if new_name != enz_type.enzyme_type:
            issues = change_enzyme_type_name(enz_type, new_name)

        return issues

    except Exception as e:
        return [f'Error updating enzyme type: {enz_type.enzyme_type}', str(e)]


def remove_enzyme_from_reactions(enz_type):
    """Remove references to enzyme type from all reactions"""

    try:
        reactions = reaction_queries.reactions_of_type([enz_type.enzyme_type])
        for reaction in reactions:
            reaction.enzyme_types.remove(enz_type.enzyme_type)
            reaction.cofactors.pop(enz_type.enzyme_type)
            reaction.save()
        return []

    except Exception as e:
        return [f'Error removing enzyme type - {enz_type.enzyme_type} from reactions', str(e)]


def delete_enzyme_type(enz_type):
    """Delete an enzyme type, if not associated with any other documents.."""

    issues = []

    seq_names = sequence_queries.seqs_of_type(enz_type.enzyme_type)
    if len(seq_names) != 0:
        for name in seq_names:
            issues.append(f'Enzyme type is present in sequence: {name}')

    reacs = reaction_queries.reactions_of_type([enz_type.enzyme_type])
    if len(reacs) != 0:
        for name in reacs:
            issues.append(f'Enzyme type is present in reaction: {name}')

    num_acts = activity_queries.activities_of_type(enz_type.enzyme_type, count_only=True)
    if num_acts != 0:
        issues.append(f"Enzyme is recorded in {num_acts} activity entries")

    if len(issues) == 0:
        try:
            remove_enzyme_from_reactions(enz_type)
            ssn_main.SSN(enz_type.enzyme_type).clear_all_data()
            enz_type.delete()
            return []
        except Exception as e:
            return ["Error deleting enzyme type", str(e)]

    else:
        return issues

def merge_enzyme_type(et_to_merge, et_merge_with):
    """Merge one enzyme type into another"""

    # 1. rename all the data
    issues = relabel_enzyme_type_data(et_to_merge, et_merge_with.enzyme_type, rename_ssn=False)
    if len(issues) != 0:
        return issues

    # 2. delete the enzyme_type being merged
    issues = delete_enzyme_type(et_to_merge)
    if len(issues) != 0:
        return issues

    return []


