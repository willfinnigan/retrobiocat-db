from pymed import PubMed

def clear_newlines_from_doi(doi):
    if doi == None:
        return 'none'

    if '\n' in doi:
        i = doi.find('\n')
        doi = doi[0:i]

    doi = doi.lower()
    doi = doi.replace(' ', '')
    return doi

def keyword_query_pubmed(keyword, date_from="1950/01/01", date_to="3000", max_results=1000):
    pubmed = PubMed(tool="RetroBioCat", email="william.finnigan@manchester.ac.uk")
    query = f'(("{date_from}"[Date - Create] : "{date_to}"[Date - Create])) AND ({keyword})'
    results = list(pubmed.query(query, max_results=max_results))
    print(f"Retrieved {len(results)} results from pubmed")

    result_list = []

    if len(results) != 0:
        for article in results:
            if article.doi != None:
                authors_list = []
                authors_string = ""
                for author_dict in article.authors:
                    authors_list.append(f"{author_dict['firstname']} {author_dict['lastname']}")
                    if authors_string != '':
                        authors_string += ', '
                    authors_string += f"{author_dict['firstname']} {author_dict['lastname']}"

                paper_dict = {'title': article.title,
                              'date': article.publication_date,
                              'journal': article.journal,
                              'short_citation': f"{article.authors[0]['lastname']} et al, {article.publication_date.year}",
                              'abstract': article.abstract,
                              'doi': clear_newlines_from_doi(article.doi),
                              'authors': authors_list}

                result_list.append(paper_dict)

    return result_list


if __name__ == '__main__':
    keyword = "carboxylic acid reductase"
    result_list = keyword_query_pubmed(keyword)

    for result in result_list:
            print(result)
