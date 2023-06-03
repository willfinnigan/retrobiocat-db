from wtforms.validators import ValidationError
import retrobiocat_web.retro.network_pathway.rdkit_utils


def is_accepted_by_rdkit(form, field):
    if retrobiocat_web.retro.network_pathway.rdkit_utils.rdkit_smile(field.data) == None:
        if field.data != '':
            raise ValidationError('SMILES not accepted by rdkit')