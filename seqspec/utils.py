import io
import gzip
from seqspec.Assay import Assay
import yaml
import requests
from Bio import GenBank


def load_spec(spec_fn: str):
    with open(spec_fn, "r") as stream:
        return load_spec_stream(stream)


def load_spec_stream(spec_stream: io.IOBase):
    data: Assay = yaml.load(spec_stream, Loader=yaml.Loader)
    # set the parent id in the Assay object upon loading it
    for r in data.library_spec:
        r.set_parent_id(None)
    return data


def load_genbank(gbk_fn: str):
    with open(gbk_fn, "r") as stream:
        return load_genbank_stream(stream)


def load_genbank_stream(gbk_stream: io.IOBase):
    data: GenBank = GenBank.read(gbk_stream)
    return data


# return cut indices for all atomic regions
def get_cuts(regions, cuts=[]):
    if not cuts:
        cuts = []
    prev = 0
    for r in regions:
        nxt = prev + r.max_len
        cuts.append((prev, nxt))
        prev = nxt
    return cuts


def write_read(header, seq, qual, f):
    f.write(f"{header}\n{seq}\n+\n{qual}\n")


def yield_onlist_contents(stream):
    for line in stream:
        yield line.strip().split()[0]


def read_list(onlist):
    """Given an onlist object read the local or remote data
    """
    if onlist.location == "remote":
        response = requests.get(onlist.filename, stream=True)
        response.raise_for_status()
        # TODO: instead of just looking at the filename, should we check the content-type?
        if onlist.filename.endswith(".gz"):
            stream = io.TextIOWrapper(gzip.GzipFile(fileobj=response.raw))
        else:
            stream = response.raw
        return list(yield_onlist_contents(stream))
    elif onlist.location == "local" and onlist.filename.endswith(".gz"):
        with gzip.open(onlist.filename, "rt") as f:
            return list(yield_onlist_contents(f))
    elif onlist.location == "local":
        with open(onlist.filename, "rt") as f:
            # just get the first column
            return list(yield_onlist_contents(f))
    else:
        raise ValueError("Unsupported location {}. Expected remote or local".format(onlist.location))


def region_ids_in_spec(seqspec, modality, region_ids):
    # return True if all region_ids are in seqspec
    spec = seqspec.get_modality(modality)
    found = []
    for region_id in region_ids:
        found += [r.region_id for r in spec.get_region_by_id(region_id)]
    return found


def file_exists(uri):
    try:
        r = requests.head(uri)
        return r.status_code == 200
    except requests.ConnectionError:
        return False


REGION_TYPE_COLORS = {
    "barcode": "#2980B9",
    "cdna": "#8E44AD",
    "fastq": "#F1C40F",
    "gdna": "#E67E22",
    "illumina_p5": "#E17A47",
    "illumina_p7": "#E17A47",
    "index5": "#4AB19D",
    "index7": "#4AB19D",
    "linker": "#1ABC9C",
    "me1": "#E74C3C",
    "me2": "#E74C3C",
    "nextera_read1": "#FF8000",
    "nextera_read2": "#FF8000",
    "poly_a": "#FF0000",
    "poly_g": "#C0C0C0",
    "poly_t": "#7F8C8D",
    "poly_c": "#2C3E50",
    "s5": "#EF3D59",
    "s7": "#EF3D59",
    "truseq_read1": "#EFC958",
    "truseq_read2": "#EFC958",
    "umi": "#16A085",
    "tag": "#344E5C",
    "protein": "#ECF0F1",
}


# unused
# '#FF8C00'
# '#95A5A6'
