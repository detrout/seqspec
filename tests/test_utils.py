import gzip
from hashlib import md5
from io import StringIO
import os
from tempfile import TemporaryDirectory
from requests import HTTPError
from unittest import TestCase
from unittest.mock import patch

from seqspec.Region import Region, Onlist

from seqspec.utils import (
    load_spec_stream, get_cuts, write_read, read_list, yield_onlist_contents
)

from .test_region import (
    region_rna_joined_dict,
    region_rna_umi_dict,
    region_rna_linker_dict,
)

example_spec = """!Assay
name: my assay
doi: https://doi.org/10.1038/nmeth.1315
publication_date: 06 April 2009
description: first method to sequence the whole transcriptome (mRNA) of a single cell
modalities:
- RNA
lib_struct: https://teichlab.github.io/scg_lib_structs/methods_html/tang2009.html
library_spec:
- !Region
  region_id: RNA
  region_type: RNA
  name: RNA
  sequence_type: joined
  sequence: CCACTACGCCTCCGCTTTCCTCTCTATGGGCAGTCGGTGATXCGCCTTGGCCGTACAGCAGNNNNNNAGAGAATGAGGAACCCGGGGCAG
  min_len: 90
  max_len: 187
  onlist: null
  regions:
  - !Region
    region_id: SOLiD_P1_adaptor
    region_type: SOLiD_P1_adaptor
    name: SOLiD_P1_adaptor
    sequence_type: fixed
    sequence: CCACTACGCCTCCGCTTTCCTCTCTATGGGCAGTCGGTGAT
    min_len: 41
    max_len: 41
    onlist: null
    regions: null
    parent_id: RNA
  - !Region
    region_id: cDNA
    region_type: cDNA
    name: cDNA
    sequence_type: random
    sequence: X
    min_len: 1
    max_len: 98
    onlist: null
    regions: null
    parent_id: RNA
  - !Region
    region_id: SOLiD_bc_adapter
    region_type: SOLiD_bc_adapter
    name: SOLiD_bc_adapter
    sequence_type: fixed
    sequence: CGCCTTGGCCGTACAGCAG
    min_len: 19
    max_len: 19
    onlist: null
    regions: null
    parent_id: RNA
  - !Region
    region_id: index
    region_type: index
    name: index
    sequence_type: onlist
    sequence: NNNNNN
    min_len: 6
    max_len: 6
    onlist: !Onlist
      filename: index_onlist.txt
      md5: null
      location: local
    regions: null
    parent_id: RNA
  - !Region
    region_id: p2_adapter
    region_type: p2_adapter
    name: p2_adapter
    sequence_type: fixed
    sequence: AGAGAATGAGGAACCCGGGGCAG
    min_len: 23
    max_len: 23
    onlist: null
    regions: null
    parent_id: RNA
"""


class TestUtils(TestCase):
    def test_load_spec_stream(self):
        with StringIO(example_spec) as instream:
            spec = load_spec_stream(instream)
        self.assertEqual(spec.name, "my assay")
        head = spec.get_modality("RNA")
        self.assertEqual(len(head.regions), 5)

    def test_get_cuts(self):
        r_umi_dict = region_rna_umi_dict("region-2")
        r_umi = Region(**r_umi_dict)
        r_linker_dict = region_rna_linker_dict("region-3")
        r_linker = Region(**r_linker_dict)
        r_expected_dict = region_rna_joined_dict("region-1", [r_umi, r_linker])
        r_expected_dict["sequence"] = r_umi_dict["sequence"] + r_linker_dict["sequence"]
        r_expected = Region(**r_expected_dict)

        r_umi_min, r_umi_max = r_umi.get_len()
        r_linker_min, r_linker_max = r_linker.get_len()
        r_linker_min += r_umi_max
        r_linker_max += r_linker_max
        cuts = get_cuts(r_expected.regions)
        self.assertEqual(cuts, [(r_umi_min, r_umi_max), (r_linker_min, r_linker_max)])

    def test_write_header(self):
        stream = StringIO()
        header = "@string"
        sequence = "CANNTG"
        quality = "IIIIII"
        write_read(header, sequence, quality, stream)

        text = stream.getvalue().split(os.linesep)
        self.assertEqual(text[0], f"{header}")
        self.assertEqual(text[1], sequence)
        self.assertEqual(text[2], "+")
        self.assertEqual(text[3], quality)

    def test_yield_onlist_contents(self):
        fake_onlist = ["ATATATAT", "GCGCGCGC"]
        fake_stream = StringIO("{}\n".format("\n".join(fake_onlist)))

        response = list(yield_onlist_contents(fake_stream))
        self.assertEqual(response, fake_onlist)

    def test_read_list_local(self):
        fake_onlist = ["ATATATAT", "GCGCGCGC"]
        fake_contents = "{}\n".format("\n".join(fake_onlist))
        fake_md5 = md5(fake_contents.encode("ascii")).hexdigest()

        with TemporaryDirectory(prefix="onlist_tmp_") as tmpdir:
            temp_list_filename = os.path.join(tmpdir, "index.txt.gz")
            with gzip.open(temp_list_filename, "wt") as stream:
                stream.write(fake_contents)

            onlist1 = Onlist(temp_list_filename, fake_md5, "local")
            loaded_list = read_list(onlist1)

            self.assertEqual(fake_onlist, loaded_list)

    def test_read_list_local_gz(self):
        fake_onlist = ["ATATATAT", "GCGCGCGC"]
        fake_contents = "{}\n".format("\n".join(fake_onlist))
        fake_md5 = md5(fake_contents.encode("ascii")).hexdigest()

        with TemporaryDirectory(prefix="onlist_tmp_") as tmpdir:
            temp_list_filename = os.path.join(tmpdir, "index.txt")
            with open(temp_list_filename, "wt") as stream:
                stream.write(fake_contents)

            onlist1 = Onlist(temp_list_filename, fake_md5, "local")
            loaded_list = read_list(onlist1)

            self.assertEqual(fake_onlist, loaded_list)

    def test_read_list_remote(self):
        fake_onlist = ["ATATATAT", "GCGCGCGC"]
        fake_contents = "{}\n".format("\n".join(fake_onlist))
        fake_md5 = md5(fake_contents.encode("ascii")).hexdigest()

        def fake_request_get(url, stream=False, **kwargs):
            class response:
                def __init__(self):
                    self.raw = StringIO(fake_contents)
                    self.status_code = 200

                def raise_for_status(self):
                    if self.status_code != 200:
                        raise HTTPError(self.status_code)

            return response()

        with patch("requests.get", new=fake_request_get):
            onlist1 = Onlist("http://localhost/testlist.txt", fake_md5, "remote")
            loaded_list = read_list(onlist1)

            self.assertEqual(fake_onlist, loaded_list)
