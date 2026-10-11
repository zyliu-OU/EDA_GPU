"""Method attribution carried with calculations, monitoring and result archives."""
ACKNOWLEDGEMENT='If EDA_GPU has benefited your research, please support us by citing this project and the GKS-EDA publications below.'
METHOD_REFERENCES=[
    dict(authors='Peifeng Su, Zhen Jiang, Zuochang Chen, Wei Wu',year=2014,
         title='Energy Decomposition Scheme Based on the Generalized Kohn–Sham Scheme',
         journal='Journal of Physical Chemistry A',volume='118',pages='2531–2542',
         doi='10.1021/jp500405s',url='https://doi.org/10.1021/jp500405s'),
    dict(authors='Peifeng Su, Zhen Tang, Wei Wu',year=2020,
         title='Generalized Kohn-Sham energy decomposition analysis and its applications',
         journal='WIREs Computational Molecular Science',article='e1460',
         doi='10.1002/wcms.1460',url='https://doi.org/10.1002/wcms.1460'),
    dict(authors='P. Su, H. Li',year=2009,
         title='Energy decomposition analysis of covalent bonds and intermolecular interactions',
         journal='Journal of Chemical Physics',volume='131',article='014102',
         doi='10.1063/1.3159673',url='https://doi.org/10.1063/1.3159673')]
TRYPTAMINE_REFERENCE=dict(authors='Zongyuan Liu, Carl O. Trindle, Quanli Gu, Wei Wu, Peifeng Su',year=2017,
    title='Unravelling hydrogen bonding interactions of tryptamine–water dimer from neutral to cation',
    journal='Physical Chemistry Chemical Physics',doi='10.1039/C7CP03491A',url='https://doi.org/10.1039/C7CP03491A')


def information(metadata=None):
    references=list(METHOD_REFERENCES)
    if metadata and 'tryptamine' in metadata.get('identity','').lower():references.append(TRYPTAMINE_REFERENCE)
    return dict(research_acknowledgement=ACKNOWLEDGEMENT,project_reference='https://github.com/zyliu-OU/EDA_GPU',method_references=references)


def text(data=None):
    data=data or information()
    return '\n'.join([data['research_acknowledgement']]+[
        f"{r['authors']} ({r['year']}). {r['title']}. DOI: {r['doi']}" for r in data['method_references']])
