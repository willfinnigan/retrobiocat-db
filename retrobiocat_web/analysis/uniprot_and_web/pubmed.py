from pymed import PubMed


def query_pubmed(doi):
    pubmed = PubMed(tool="RetroBioCat", email="william.finnigan@manchester.ac.uk")
    query = f'{doi}[Location ID]'
    results = list(pubmed.query(query, max_results=1))

    title = ''
    authors_list = []
    date = ''
    journal = ''
    cite_mini = ''

    if (len(results) != 0):
        article = results[0]
        if article.doi == doi:
            title = article.title
            journal = article.journal
            authors_list = []
            for author_dict in article.authors:
                authors_list.append(f"{author_dict['firstname']} {author_dict['lastname']}")
            date = article.publication_date
            cite_mini = f"{article.authors[0]['lastname']} et al, {date.year}"
    return title, authors_list, journal, date, cite_mini