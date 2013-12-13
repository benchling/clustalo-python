from multiprocessing import Pipe, Process

from _clustalo import DNA, RNA, PROTEIN, clustalo as clustalo_c

__all__ = ['DNA', 'RNA', 'PROTEIN', 'clustalo']

def _clustalo_helper(conn, data, **kwargs):
    results = clustalo_c(data, **kwargs)
    conn.send(results)

def clustalo(data, timeout=None, seqtype=DNA,
             mbed_guide_tree=True,
             mbed_iteration=True,
             num_combined_iterations=0,
             max_guidetree_iterations=-1,
             max_hmm_iterations=-1,
             num_threads=1):
    """
    Runs clustal omega in a separate process.
    Returns a dictionary of sequence_named => aligned_bases ('_' for gaps)
    If timeout is reached, returns None

    Arguments:
    data -- dictionary of sequence_name => bases

    Optional arguments:
    timeout -- maximum time (seconds) before aborting
    seqtype -- should be one of clustalo.DNA, clustalo.RNA, or clustalo.PROTEIN
    mbed_guide_tree -- whether mBed-like clustering guide tree should be used
    mbed_iteration -- whether mBed-like clustering iteration should be used
    num_combined_iterations -- number of (combined guide-tree/HMM) iterations
    max_guidetree_iterations -- max guide tree iterations within combined iterations
    max_hmm_iterations -- max HMM iterations within combined iterations
    num_threads -- number of threads to use (requires libclustalo compiled with OpenMP)
    """
    kwargs = {
        'seqtype': seqtype,
        'mbed_guide_tree': mbed_guide_tree,
        'mbed_iteration': mbed_iteration,
        'num_combined_iterations': num_combined_iterations,
        'max_guidetree_iterations': max_guidetree_iterations,
        'max_hmm_iterations': max_hmm_iterations,
        'num_threads': num_threads
    }

    parent_conn, child_conn = Pipe()
    p = Process(target=_clustalo_helper, args=(child_conn, data,), kwargs=kwargs)
    results = None
    try:
        p.start()
        if parent_conn.poll(timeout):
            results = parent_conn.recv()
            p.join()
    finally:
        p.terminate()
        p.join()

    return results
