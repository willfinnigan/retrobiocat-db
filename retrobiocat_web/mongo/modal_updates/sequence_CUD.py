from retrobiocat_web.mongo.models.biocatdb_models import Sequence, EnzymeType, UniRef50, Activity, Tag, SequenceOtherNamesData
from retrobiocat_web.mongo.functions.save_sequence_data import save_sequence_functions
from distutils.util import strtobool
from retrobiocat_web.mongo.modal_updates import bioinformatics_CUD
from retrobiocat_web.mongo.model_queries import activity_queries, paper_queries, sequence_queries
import json

"""
Use only these functions for creating, updating or deleting Sequence objects.
Each function will return list of issues (which is empty if everything worked fine)
These functions keep related fields and bioinformatics files up to date.
"""

def create_new_sequence(seq_dict, user=None, paper=None):
    """
    Function to create a new sequence
    Returns both issues and the seq_obj is successful
    """

    issues = save_sequence_functions.check_seq_dict(seq_dict)
    if len(issues) != 0:
        return issues, None

    try:
        sequence = save_sequence_functions.sanitise_string(seq_dict.get('sequence', ''))
        sequence = save_sequence_functions.if_sequence_is_DNA_then_convert_to_protein(sequence)
        name = save_sequence_functions.sanitise_string(seq_dict.get('enzyme_name', ''))
        other_names = []
        if seq_dict.get('other_names', '') != '':
            other_names = seq_dict.get('other_names', '').replace(' ', '').split(', ')
            for name in other_names:
                if save_sequence_functions.sequence_name_check(name) == True:
                    return ['Sequence with this name already exists (either as the name, or in other names)'], None

        if save_sequence_functions.sequence_name_check(name) == True:
            return ['Sequence with this name already exists (either as the name, or in other names)'], None

        for name in other_names:
            if save_sequence_functions.sequence_name_check(name) == True:
                return ['Sequence with this name already exists (either as the name, or in other names)'], None

        # force sequence_unavailable into the correct format
        seq_unavailable = str(seq_dict.get('sequence_unavailable', 'false')).lower()
        if seq_unavailable in ['1', '1.0', 'yes', 'true']:
            seq_dict['sequence_unavailable'] = True
        elif seq_unavailable in ['0', '0.0', 'no', 'false']:
            seq_dict['sequence_unavailable'] = False

        papers = []
        if paper is not None:
            papers = [paper]

        seq = Sequence(enzyme_name=name,
                       enzyme_type=seq_dict['enzyme_type'],
                       other_names_data=[SequenceOtherNamesData(name=name) for name in other_names],
                       sequence=sequence,
                       n_tag=seq_dict.get('n_tag', ''),
                       c_tag=seq_dict.get('c_tag', ''),
                       sequence_unavailable=bool(strtobool(str(seq_dict.get('sequence_unavailable', 'False')))),
                       accession=seq_dict.get('accession', ''),
                       other_identifiers=seq_dict.get('other_names', '').split(', '),
                       pdb=seq_dict.get('pdb', ''),
                       mutant_of=seq_dict.get('mutant_of', ''),
                       notes=seq_dict.get('notes', ''),
                       bioinformatics_ignore=bool(strtobool(str(seq_dict.get('bioinformatics_ignore', 'False')))),
                       owner=user,
                       added_by=user,
                       papers=papers)
        seq.save()

        issues += find_and_update_sequence_tags(seq)

        # return issues and sequence object
        return issues, seq

    except Exception as e:
        # return the exception as entry in list of issues
        return [str(e)], None

def update_sequence(seq_obj, new_sequence_string):
    """Update the sequence field"""

    seq_string = save_sequence_functions.sanitise_sequence(new_sequence_string)
    bad_chars = save_sequence_functions.sequence_check(seq_string)
    if len(bad_chars) != 0:
        return [f'Invalid sequence chars: {bad_chars}']

    seq_obj.sequence = seq_string
    seq_obj.save()
    issues = find_and_update_sequence_tags(seq_obj)

    seq_needs_new_embedding(seq_obj)
    seq_needs_realigning_in_ssn(seq_obj)
    seq_needs_re_blasting(seq_obj)

    return issues

def update_enzyme_name(seq_obj, new_name):
    """Update the enzyme_name field"""
    new_name = save_sequence_functions.sanitise_string(new_name)
    old_name = seq_obj.enzyme_name
    if old_name == new_name:
        return ['Name already exists']

    if new_name in seq_obj.other_names:
        seq_obj.other_names.remove(new_name)
        seq_obj.save()

    if save_sequence_functions.sequence_name_check(new_name) == True:
        return ['Sequence with this name already exists (either as the name, or in other names)'], None

    print(f"Updating sequence name: {new_name}")

    # update mutants_of
    mutants = Sequence.objects(mutant_of=seq_obj.enzyme_name)
    for mut in mutants:
        mut.mutant_of = new_name
        mut.save()

    # update activity
    acts = Activity.objects(enzyme_name=seq_obj.enzyme_name)
    for act in acts:
        act.enzyme_name = new_name
        act.save()

    # update to the new name
    seq_obj.enzyme_name = new_name
    seq_obj.save()

    # if name is changes it'll have to be re-added to the ssn
    seq_needs_realigning_in_ssn(seq_obj)

    return []

def update_enzyme_type(seq_obj, new_type):
    """Update the enzyme_type field"""

    old_type = seq_obj.enzyme_type

    # new_type must exist in EnzymeType collection
    if EnzymeType.objects(enzyme_type=new_type).count() == 0:
        return ['Enzyme type not found']

    # don't bother updating if no change
    if old_type == new_type:
        return ['Enzyme is already of this type']

    # mark the bioinformatics for the old_type as needing updating
    seq_needs_realigning_in_ssn(seq_obj)
    seq_needs_re_blasting(seq_obj)

    # update activity
    acts = Activity.objects(enzyme_name=seq_obj.enzyme_name)
    for act in acts:
        act.enzyme_type = new_type
        act.save()

    # change the type and save
    seq_obj.enzyme_type = new_type
    seq_obj.save()

    # mark bioinformatics of new type as needing updating also
    seq_needs_realigning_in_ssn(seq_obj)
    seq_needs_re_blasting(seq_obj)

    return []

def create_other_name(seq_obj, other_names_dict):
    new_name = save_sequence_functions.sanitise_string(other_names_dict.get('new_name', None))
    if new_name is None or new_name == '':
        return [f"Must enter a name"]
    if save_sequence_functions.sequence_name_check(new_name) == True:
       return [f"Can create other name {new_name} as this name already exists (either as an other name or a sequence name)"]

    try:
        other_name_data = SequenceOtherNamesData(name=new_name,
                                                 n_tag=other_names_dict.get('n_tag', None),
                                                 c_tag=other_names_dict.get('c_tag', None),
                                                 notes=other_names_dict.get('notes', None))
        seq_obj.other_names_data.append(other_name_data)
        seq_obj.save()
    except Exception as e:
        return [f'Error creating other name {new_name}', str(e)]

    return []

def update_other_name(seq_obj, other_names_dict):
    """ Update a sequences other names data field, ensure no duplicates"""

    issues = []

    existing_name = other_names_dict.get('existing_name', None)
    if existing_name is None or existing_name == '':
        return ['Error - no name provided for updating other names']

    # loop over current other names, if exists update
    for i, data in enumerate(seq_obj.other_names_data):
        if data.name == existing_name:
            print('name matches')
            # update name
            new_name = other_names_dict.get('new_name', None)
            if new_name is not None and new_name != existing_name:
                new_name = save_sequence_functions.sanitise_string(new_name)
                if save_sequence_functions.sequence_name_check(new_name) == True:
                    return [f"Can not update name to {new_name} as this name already exists (either as an other name or a sequence name)"]
                else:
                    seq_obj.other_names_data[i].name = new_name

            # update all other fields
            try:
                seq_obj.other_names_data[i].n_tag = other_names_dict.get('n_tag', None)
                seq_obj.other_names_data[i].c_tag = other_names_dict.get('c_tag', None)
                seq_obj.other_names_data[i].notes = other_names_dict.get('notes', None)
                seq_obj.save()
                return issues
            except Exception as e:
                return [
                    f'Error updating other_names for sequence {seq_obj.enzyme_name} - other name {existing_name}',
                    str(e)]

    return [f'Error - {existing_name} did not match a record, nothing to update']








def delete_other_name(seq_obj, other_name):
    ''' Delete a other name entry matching the name'''

    # loop over current other names, if exists delete, otherwise error
    for i, data in enumerate(seq_obj.other_names_data):
        if data.name == other_name:
            seq_obj.other_names_data.pop(i)
            seq_obj.save()
            return []

    return ['No other name matches, nothing deleted']

def update(seq_obj, seq_dict, user=None):
    """Update a sequence entry from the given seq_dict"""

    issues = []

    # if the sequence is marked as reviewed, then cannot use this function
    if seq_obj.reviewed == True:
        issues = ['Sequence is marked as reviewed, cannot update essential fields in this way']
        return issues

    # update metadata fields
    try:

        if 'sequence_unavailable' in seq_dict:
            seq_obj.sequence_unavailable = bool(strtobool(str(seq_dict['sequence_unavailable'])))
        if 'other_identifiers' in seq_dict:
            seq_obj.other_identifiers = seq_dict['other_identifiers'].split(', ')
        if 'bioinformatics_ignore' in seq_dict:
            seq_obj.bioinformatics_ignore = bool(strtobool(str(seq_dict['bioinformatics_ignore'])))
        if 'n_tag' in seq_dict:
            seq_obj.n_tag = seq_dict['n_tag']
        if 'c_tag' in seq_dict:
            seq_obj.c_tag = seq_dict['c_tag']
        if 'accession' in seq_dict:
            seq_obj.accession = seq_dict['accession']
        if 'pdb' in seq_dict:
            seq_obj.pdb = seq_dict['pdb']
        if 'notes' in seq_dict:
            seq_obj.notes = seq_dict['notes']
        if 'mutant_of' in seq_dict:
            seq_obj.mutant_of = seq_dict['mutant_of']

        save_sequence_functions.clean_seq_lists(seq_obj)

        seq_obj.save()

        if user is not None:
            add_edited_by(seq_obj, user)

    except Exception as e:
        return ['Failed to update other sequence metadata', str(e)]

    # if there's a new enzyme_type, update this
    if seq_dict.get('enzyme_type', '') not in ['', seq_obj.enzyme_type]:
        issues = update_enzyme_type(seq_obj, seq_dict['enzyme_type'])
        if len(issues) != 0:
            return issues

    # if there's a new name, update this
    if seq_dict.get('enzyme_name', '') not in ['', seq_obj.enzyme_name]:
        issues = update_enzyme_name(seq_obj, seq_dict['enzyme_name'])
        if len(issues) != 0:
            return issues

    # if there's a new sequence, update this
    if seq_dict.get('sequence', '') not in ['', seq_obj.sequence]:
        issues = update_sequence(seq_obj, seq_dict['sequence'])
        if len(issues) != 0:
            return issues

    # if all went smoothly issues is a blank list
    return issues

def update_non_essential_only(seq_obj, seq_dict, user=None):
    """If a sequence has been reviewed, key fields are locked so use this function to update instead"""

    issues = []

    # update other fields
    try:
        if 'other_identifiers' in seq_dict:
            seq_obj.other_identifiers = seq_dict['other_identifiers'].split(', ')
        if 'pdb' in seq_dict:
            seq_obj.pdb = seq_dict['pdb']
        if 'notes' in seq_dict:
            seq_obj.notes = seq_dict['notes']

        save_sequence_functions.clean_seq_lists(seq_obj)

        seq_obj.save()

        if user is not None:
            add_edited_by(seq_obj, user)

        return issues

    except Exception as e:
        return ['Failed to update non-essential sequence metadata', str(e)]

def change_owner(seq_obj, user):
    """ Update the owner of a sequence"""

    try:
        seq_obj.owner = user
        seq_obj.save()
        return []
    except Exception as e:
        return [f"Couldn't update user for sequence: {seq_obj.enzyme_name}", str(e)]

def add_edited_by(seq_obj, user):
    if user is None:
        return

    if user not in seq_obj.edits_by:
        seq_obj.edits_by.append(user)

def update_reviewed(seq_obj, review_status, check=True, update_bio=True):
    """Update reviewed.  If reviewed is true need to check first"""
    issues = []
    if review_status == True and check == True:
        if not save_sequence_functions.check_seq_has_protein_seq_or_unavailable(seq_obj):
            return ['Sequence must have a protein sequence or be marked as sequence unavailable before reviewing.']

    try:
        seq_obj.reviewed = review_status
        seq_obj.save()
    except Exception as e:
        return [f'Error updating reviewed status for {seq_obj.enzyme_name}', str(e)]

    if update_bio == True:
        bioinformatics_CUD.enzyme_type_needs_blasting(seq_obj.enzyme_type)
        bioinformatics_CUD.SSN_needs_updating(seq_obj.enzyme_type)

    return issues

def add_paper(seq_obj, paper, allow_while_reviewed=False):
    """Add paper to sequence entry"""

    # if paper.seq_reviewed is true, cannot add any papers
    if paper.seq_reviewed == True and allow_while_reviewed == False:
        return ['The sequences section of this paper has been reviewed, no new sequences can be added']

    try:
        if paper not in seq_obj.papers:
            seq_obj.papers.append(paper)
            seq_obj.save()
            return []
        else:
            return ['Sequence has already been added to this paper']
    except Exception as e:
        return [f'Error adding papers to {seq_obj.enzyme_name}', str(e)]

def remove_paper(seq_obj, paper):
    """Remove paper from sequence entry"""

    # if paper.seq_reviewed is true, cannot remove any papers
    if paper.seq_reviewed == True:
        return ['The sequences section of this paper has been reviewed, no new sequences can be removed']

    if paper not in seq_obj.papers:
        return ['Sequence is not in paper']

    num_act = activity_queries.activity_in_paper(paper, seq_name=seq_obj.enzyme_name, count_only=True)
    if num_act != 0:
        return ['This sequence still has activity records in this paper, remove first']

    try:
        seq_obj.papers.remove(paper)
        seq_obj.save()

        issues = []
        if len(seq_obj.papers) == 0:
            issues = delete_sequence(seq_obj)

        return issues

    except Exception as e:
        return [f'Error removing papers to {seq_obj.enzyme_name}', str(e)]

def seq_needs_re_blasting(seq_obj):
    # remove this sequence as a source for uniref50 hits
    uniref_query = UniRef50.objects(result_of_blasts_for=seq_obj)
    for uniref in uniref_query:
        uniref.result_of_blasts_for.remove(seq_obj)
        uniref.save()
    seq_obj.blast = None
    bioinformatics_CUD.enzyme_type_needs_blasting(seq_obj.enzyme_type)

def seq_needs_realigning_in_ssn(seq_obj):
    # Mark sequence as needing alignments made on ssn
    seq_obj.alignments_made = False
    seq_obj.save()
    bioinformatics_CUD.SSN_needs_updating(seq_obj.enzyme_type)

def seq_needs_new_embedding(seq_obj):
    seq_obj.seqvec = None
    seq_obj.unirep = None
    seq_obj.save()

def update_unirep(seq_obj, embedding):
    seq_obj.unirep = json.dumps(list(embedding))
    seq_obj.save()

def delete_sequence(seq_obj):
    """Delete a sequence, but only if no associated data"""

    issues = []

    # Do not delete if there is activity data
    if Activity.objects(enzyme_name=seq_obj.enzyme_name).count() != 0:
        issues.append(f"Sequence is recorded in activity data for papers...")
        papers = []
        for act in Activity.objects(enzyme_name=seq_obj.enzyme_name):
            if act.short_citation not in papers:
                papers.append(act.short_citation)
        for paper in papers:
            issues.append(f"Sequence is recorded in activity data for {paper}")

    # Do not delete if there are mutants of this enzyme
    if Sequence.objects(mutant_of=seq_obj.enzyme_name).count() != 0:
        issues.append('There are mutants of this enzymes...')
        for mut in Sequence.objects(mutant_of=seq_obj.enzyme_name):
            issues.append(f"Sequence is a parent of mutant {mut.enzyme_name}")

    # Delete if no issues
    if len(issues) == 0:
        seq_obj.delete()
        return []

    return issues

def merge_sequences(seq_obj_merge_to, seq_obj_merge_from, save_other_name=True):
    issues = []

    # check merging is ok
    if seq_obj_merge_to.enzyme_name == seq_obj_merge_from.enzyme_name:
        issues.append('Enzymes are the same, cannot merge')
    if seq_obj_merge_to.enzyme_type != seq_obj_merge_from.enzyme_type:
        issues.append('Enzymes are not the same type, cannot merge')
    if len(issues) != 0:
        return issues

    # -- perform merge --
    try:

        # add papers
        for paper in seq_obj_merge_from.papers:
            add_paper(seq_obj_merge_to, paper, allow_while_reviewed=True)

        # rename activity
        acts = Activity.objects(enzyme_name=seq_obj_merge_from.enzyme_name)
        for act in acts:
            act.enzyme_name = seq_obj_merge_to.enzyme_name
            act.save()

        # rename mutant_of fields
        mutants = Sequence.objects(mutant_of=seq_obj_merge_from.enzyme_name)
        for mut in mutants:
            mut.mutant_of = seq_obj_merge_to.enzyme_name
            mut.save()

        # only copy other name data over if true
        if save_other_name == True:

            # copy over any existing other names documents:
            seq_obj_merge_to.other_names_data += seq_obj_merge_from.other_names_data

            # add old name to other names
            papers = []
            if seq_obj_merge_from.papers is not None:
                papers = [str(paper.id) for paper in seq_obj_merge_from.papers]
            new_other_name_dict = {'existing_name': seq_obj_merge_from.enzyme_name,
                                   'n_tag': seq_obj_merge_from.n_tag,
                                   'c_tag': seq_obj_merge_from.c_tag,
                                   'notes': seq_obj_merge_from.c_tag,
                                   'papers': papers}

            # remove any blank values
            new_other_name_dict = {k: v for k, v in new_other_name_dict.items() if v != ''}

            # create new other_names_data documents
            issues += update_other_names(seq_obj_merge_to, new_other_name_dict)

        # final save and delete
        seq_obj_merge_to.save()
        seq_obj_merge_from.delete()

        return issues

    except Exception as e:
        return ['Error merging sequences', str(e)]



def find_and_update_sequence_tags(seq_obj, n_tags=None, c_tags=None):
    """
    If a sequence has one of the predefined sequence tags at the N or C term,
    then remove this and add it to the n_tag or c_tag field
    """

    # if theres no protein sequence, nothing to do..
    if seq_obj.sequence == "" or seq_obj.sequence == None:
        return []

    try:
        # load the tags if not given
        # important to sort in terms of length so longer tags are tested first..
        if n_tags == None:
            n_tags = Tag.objects(n_term=True).distinct('seq')
            n_tags = sorted(n_tags, key=len, reverse=True)
        if c_tags == None:
            c_tags = Tag.objects(c_term=True).distinct('seq')
            c_tags = sorted(c_tags, key=len, reverse=True)

        # does sequence have a tag? if so remove and add to appropriate tag field
        for n_tag in n_tags:
            if n_tag == seq_obj.sequence[0:len(n_tag)]:
                print(n_tag)
                seq_obj.n_tag = n_tag
                print(seq_obj.n_tag)
                seq_obj.sequence = seq_obj.sequence[len(n_tag):]
                seq_obj.save()
        for c_tag in c_tags:
            if c_tag == seq_obj.sequence[0:len(c_tag)]:
                seq_obj.c_tag = c_tag
                seq_obj.sequence = seq_obj.sequence[len(c_tag):]
                seq_obj.save()

        # return no issues
        return []

    except Exception as e:
        return [f'Error while trying to identify tags in protein sequence for {seq_obj.enzyme_name}', str(e)]


def update_other_names_paper_use(alt_names, paper_id):
    try:
        paper = paper_queries.paper_from_id(paper_id)
        paper_seqs = sequence_queries.seqs_of_paper(paper)
        for seq in paper_seqs:
            for i, other_names_data in enumerate(seq.other_names_data):
                # first remove ref to paper in other_names data
                if paper in other_names_data.papers:
                    seq.other_names_data[i].papers.remove(paper)

                # then set add back if alt name is given
                if other_names_data.name in alt_names:
                    seq.other_names_data[i].papers.append(paper)

                seq.save()
        return []

    except Exception as e:
        return [f"Error updating other names", str(e)]


if __name__ == '__main__':
    seq_dict = {}
    test = bool(strtobool(str(seq_dict.get('sequence_unavailable', 'False'))))
    print(test)


    from retrobiocat_web.mongo.default_connection import make_default_connection
    make_default_connection()

    seq_obj = Sequence.objects(enzyme_name='BsLeuDH').first()


    print(seq_obj.to_json())

    other_names = Sequence.objects(enzyme_name='mpCAR').distinct('other_names_data.name')
    print(other_names)



