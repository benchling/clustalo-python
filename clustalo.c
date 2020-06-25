#include <Python.h>
#include <clustal-omega.h>

static PyObject *
clustalo_clustalo(PyObject *self, PyObject *args, PyObject *keywds)
{
    mseq_t *prMSeq = NULL;

    LogDefaultSetup(&rLog);
    LogMuteAll(&rLog);

    // Initialize sequence and alignment options.
    NewMSeq(&prMSeq);
    // Assuming input sequences are not aligned.
    prMSeq->aligned = false;

    opts_t rAlnOpts;
    SetDefaultAlnOpts(&rAlnOpts);

    // Required
    PyObject *inputDict;
    // Optional
    int seqtype = SEQTYPE_DNA;
    PyObject *mbedGuideTree = NULL;
    PyObject *mbedIteration = NULL;
    int numCombinedIterations = rAlnOpts.iNumIterations;
    int maxGuidetreeIterations = rAlnOpts.iMaxGuidetreeIterations;
    int maxHMMIterations = rAlnOpts.iMaxHMMIterations;
    int numThreads = 1;
    int outOrder = 1;
    static char *kwlist[] = {
        "seqs",
        "seqtype",
        "mbed_guide_tree",
        "mbed_iteration",
        "num_combined_iterations",
        "max_guidetree_iterations",
        "max_hmm_iterations",
        "num_threads",
        "output_order",
        NULL
    };
    if (!PyArg_ParseTupleAndKeywords(args, keywds, "O!|iOOiiiii", kwlist,
            &PyDict_Type, &inputDict,
            &seqtype,
            &mbedGuideTree,
            &mbedIteration,
            &numCombinedIterations,
            &maxGuidetreeIterations,
            &maxHMMIterations,
            &numThreads,
            &outOrder))
        return NULL;

    if (PyObject_Not(inputDict))
        return PyDict_New();

    InitClustalOmega(numThreads);

    switch (seqtype) {
        case SEQTYPE_DNA:
        case SEQTYPE_RNA:
        case SEQTYPE_PROTEIN:
            prMSeq->seqtype = seqtype;
            break;
        default:
            PyErr_SetString(PyExc_ValueError, "Bad seqtype, must be one of SEQTYPE_DNA, SEQTYPE_RNA, SEQTYPE_PROTEIN, or SEQTYPE_UNKNOWN");
            return NULL;
    }
    if (mbedGuideTree != NULL)
        rAlnOpts.bUseMbed = PyObject_IsTrue(mbedGuideTree);
    if (mbedIteration != NULL)
        rAlnOpts.bUseMbedForIteration = PyObject_IsTrue(mbedIteration);
    rAlnOpts.iNumIterations = numCombinedIterations;
    rAlnOpts.iMaxGuidetreeIterations = maxGuidetreeIterations;
    rAlnOpts.iMaxHMMIterations = maxHMMIterations;

    // allocating the tree_order of prMSeq is enough to capture the tree_order information
    if (outOrder == 1) {
        prMSeq->tree_order = (int *)CKMALLOC(prMSeq->nseqs * sizeof(int));
    }
    rAlnOpts.iOutputOrder = outOrder;

    // Read in sequences from input.
    PyObject *key, *value;
    Py_ssize_t pos = 0;
    while (PyDict_Next(inputDict, &pos, &key, &value)) {
        #if PY_MAJOR_VERSION >= 3
            char *seq = PyBytes_AsString(PyUnicode_AsUTF8String(value));
        #else
            char *seq = PyString_AsString(value);
        #endif
        // Sanitize sequence.
        int seqPos;
        for (seqPos = 0; seqPos < (int)strlen(seq); seqPos++) {
            char *res = &(seq[seqPos]);
            if (isgap(*res))
                continue;
            if (prMSeq->seqtype == SEQTYPE_PROTEIN && strchr(AMINO_ALPHABET, toupper(*res)) == NULL) {
                *res = AMINOACID_ANY;
            } else if (prMSeq->seqtype==SEQTYPE_DNA && strchr(DNA_ALPHABET, toupper(*res)) == NULL) {
                *res = NUCLEOTIDE_ANY;
            } else if (prMSeq->seqtype==SEQTYPE_RNA && strchr(RNA_ALPHABET, toupper(*res)) == NULL) {
                *res = NUCLEOTIDE_ANY;
            }
        }
        #if PY_MAJOR_VERSION >= 3
            AddSeq(&prMSeq, PyBytes_AsString(PyUnicode_AsUTF8String(key)), seq);
        #else
            AddSeq(&prMSeq, PyString_AsString(key), seq);
        #endif
    }
    // Can't align with only 1 sequence.
    if (prMSeq->nseqs <= 1) {
        return PyDict_Copy(inputDict);
    }

    // Perform the alignment.
    int rv;
    Py_BEGIN_ALLOW_THREADS
    rv = Align(prMSeq, NULL, &rAlnOpts);
    Py_END_ALLOW_THREADS

    if (rv) {
        PyErr_SetString(PyExc_RuntimeError, "Error while running clustal omega alignment.");
        return NULL;
    }

    // Return the aligned results in a dict.
    PyObject *returnDict = PyDict_New();
    int idx;
    if (outOrder == 1){
        for (idx = 0; idx < prMSeq->nseqs; idx++) {
            //printf("NAME OF SEQUENCE: %s, %i \n", prMSeq->sqinfo[prMSeq->tree_order[idx]].name, prMSeq->tree_order[idx]);
            const char *key = prMSeq->sqinfo[prMSeq->tree_order[idx]].name;
#if PY_MAJOR_VERSION >= 3
            PyObject *value = PyUnicode_FromString(prMSeq->seq[prMSeq->tree_order[idx]]);
#else
            PyObject *value = PyString_FromString(prMSeq->seq[prMSeq->tree_order[idx]]);
#endif
            PyDict_SetItemString(returnDict, key, value);
        }
    }
    else {
        for (idx = 0; idx < prMSeq->nseqs; idx++) {
            const char *key = prMSeq->sqinfo[idx].name;
#if PY_MAJOR_VERSION >= 3
            PyObject *value = PyUnicode_FromString(prMSeq->seq[idx]);
#else
            PyObject *value = PyString_FromString(prMSeq->seq[idx]);
#endif
            PyDict_SetItemString(returnDict, key, value);
        }
    }
    return returnDict;
}
#if PY_MAJOR_VERSION >= 3
  #define MOD_ERROR_VAL NULL
  #define MOD_SUCCESS_VAL(val) val
  #define MOD_INIT(name) PyMODINIT_FUNC PyInit_##name(void)
  #define MOD_DEF(ob, name, doc, methods) \
          static struct PyModuleDef moduledef = { \
            PyModuleDef_HEAD_INIT, name, doc, -1, methods, }; \
          ob = PyModule_Create(&moduledef);
#else
  #define MOD_ERROR_VAL
  #define MOD_SUCCESS_VAL(val)
  #define MOD_INIT(name) void init##name(void)
  #define MOD_DEF(ob, name, doc, methods) \
          ob = Py_InitModule3(name, methods, doc);
#endif


static PyMethodDef ClustaloMethods[] = {
    {"clustalo",  (PyCFunction)clustalo_clustalo, METH_VARARGS | METH_KEYWORDS,
     "Runs clustal omega.\n"
     "\n"
     "Args:\n"
     "  data (dict): dictionary of sequence_name => bases\n"
     "\n"
     "Kwargs:\n"
     "  seqtype (int): should be one of clustalo.DNA, clustalo.RNA, or clustalo.PROTEIN\n"
     "  mbed_guide_tree (bool): whether mBed-like clustering guide tree should be used\n"
     "  mbed_iteration (bool): whether mBed-like clustering iteration should be used\n"
     "  num_combined_iterations (int): number of (combined guide-tree/HMM) iterations\n"
     "  max_guidetree_iterations (int): max guide tree iterations within combined iterations\n"
     "  max_hmm_iterations (int): max HMM iterations within combined iterations\n"
     "  num_threads (int): number of threads to use (requires libclustalo compiled with OpenMP)\n"
     "  output_order (int): return the alignment with either the input order (0) or alignment tree order (1)\n"
     "\n"
     "Returns dict of sequence_names:aligned_bases ('-' for gaps)\n"},
    {NULL, NULL, 0, NULL}
};

MOD_INIT(clustalo)
{
    PyObject *module;
    MOD_DEF(module, "clustalo", NULL, ClustaloMethods);
    if (module == NULL) {
        return MOD_ERROR_VAL;
    }
    PyModule_AddIntConstant(module, "DNA", SEQTYPE_DNA);
    PyModule_AddIntConstant(module, "RNA", SEQTYPE_RNA);
    PyModule_AddIntConstant(module, "PROTEIN", SEQTYPE_PROTEIN);

    return MOD_SUCCESS_VAL(module);
}
