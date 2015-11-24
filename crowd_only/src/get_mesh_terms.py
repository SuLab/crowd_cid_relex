# Tong Shu Li
# Last updated: 2015-11-24
"""Grab MeSH terms for a PubMed abstract using Eutils."""

from collections import defaultdict
import numpy as np
import pandas as pd
import requests
import xml.etree.cElementTree as ET

def query_ncbi(eutil, params):
    BASE = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
    resp = requests.get(BASE + eutil, params = params)

    assert resp.status_code == requests.codes.ok, "Response code #{}".format(
        resp.status_code
    )
    return resp.text

def get_pubmed_article_xml_tree(pubmed_id):
    params = {
        "db": "pubmed",
        "id": pubmed_id,
        "rettype": "abstract"
    }

    resp = query_ncbi("efetch.fcgi", params)
    return ET.fromstring(resp)


class Article:
    """A PubMed article summary.

    Given a PubMed identifier, queries PubMed's API to grab relevant info.
    """
    def __init__(self, pmid):
        """Creates an article summary by querying the PubMed API."""

        def get_publication_type():
            """Parse out publication type classifications"""
            temp = defaultdict(list)
            for node in tree.findall(".//PublicationTypeList/PublicationType"):
                temp["pub_type"].append(node.attrib["UI"])
                temp["pub_name"].append(node.text)

            return pd.DataFrame(temp)

        def get_mesh_terms():
            """Parse out MeSH annotations from the XML and store as a dataframe.
            """
            temp = defaultdict(list)
            for node in tree.findall(".//MeshHeadingList/MeshHeading"):
                desc = node.find("DescriptorName")
                desc_major = desc.attrib["MajorTopicYN"] == "Y"

                qualifiers = node.findall("QualifierName")
                if not qualifiers:
                    temp["descriptor"].append(desc.attrib["UI"])
                    temp["desc_name"].append(desc.text)
                    temp["qualifier"].append(np.nan)
                    temp["qual_name"].append(np.nan)
                    temp["major_topic"].append(desc_major)
                else:
                    for qual in qualifiers:
                        temp["descriptor"].append(desc.attrib["UI"])
                        temp["desc_name"].append(desc.text)
                        temp["qualifier"].append(qual.attrib["UI"])
                        temp["qual_name"].append(qual.text)
                        temp["major_topic"].append(desc_major or qual.attrib["MajorTopicYN"] == "Y")

            return pd.DataFrame(temp)

        self.pmid = int(pmid)

        tree = get_pubmed_article_xml_tree(pmid)
        loc = tree.find("./PubmedArticle/MedlineCitation/Article")

        # get the title
        node = loc.find("./ArticleTitle")
        assert node is not None, "No title for PMID:{}!".format(pmid)
        self.title = node.text

        # get the abstract (which may not exist)
        node = loc.find("./Abstract/AbstractText")
        self.abstract = node.text if node is not None else None

        # parse out all the chemicals mentioned
        self.chemical_terms = [node.attrib["UI"]
            for node in tree.findall(".//ChemicalList/Chemical/NameOfSubstance")
        ]

        # parse out the mesh annotations
        self.mesh_terms = get_mesh_terms()

        # get the type of publication
        self.pub_type = get_publication_type()

    def __repr__(self):
        return "<{}>: PMID:{}. '{}'\n{} MeSH terms".format(self.__class__.__name__,
            self.pmid, self.title, len(self.mesh_terms))
