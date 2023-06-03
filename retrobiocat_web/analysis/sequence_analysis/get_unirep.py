from retrobiocat_web.mongo.modal_updates import sequence_CUD
from retrobiocat_web.mongo.model_queries import sequence_queries
from jax_unirep import get_reps
# https://elarkk.github.io/jax-unirep/getting-started/

def remove_stop_codons(seq):
    return seq.replace('*', '')


def generate_embedding(seq):
    try:
        protein_seq = str(seq.replace("*", '').replace(' ',''))
        h_avg, h_final, c_final = get_reps(protein_seq)
        embedding = h_avg.tolist()
        return embedding[0]
    except Exception as e:
        print(e)
        return None

def generate_and_save_embedding(list_seq_names, log=False):
    """For a list of seq names, generate the embeddings and save"""

    seq_objs = sequence_queries.seq_objs_from_seq_names(list_seq_names)

    # generate embeddings and save
    for seq_obj in seq_objs:
        emb = generate_embedding(seq_obj.sequence)
        if log == True:
            print(emb)
        if emb is None:
            print(f'Error generating embedding for {seq_obj.enzyme_name}')
            print(f"Sequence = {seq_obj.sequence}")
        else:
            sequence_CUD.update_unirep(seq_obj, emb)


if __name__ == '__main__':

    import mongoengine as db
    from retrobiocat_web.mongo.default_connection import make_default_connection
    from retrobiocat_web.mongo.models.biocatdb_models import Sequence

    make_default_connection()

    seqs = Sequence.objects(db.Q(enzyme_type='CAR') & db.Q(sequence__ne=""))

    generate_and_save_embedding(['mpCAR'], log=True)

