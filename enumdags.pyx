"""Main code for classifying GDAGS"""

# see Joe Henson, Raymond Lal, Matthew F. Pusey
# "Theory-independent limits on correlations from generalised Bayesian
#  networks", arXiv:1405.2572

# Copyright (C) 2014 Matthew F. Pusey -- see LICENSE.txt

import itertools, multiprocessing
from bitfuncs import *
from closure import *

from libc.stdlib cimport malloc, free
from libc.string cimport memset

cdef extern from *:
    int _ctz "__builtin_ctzll" (unsigned long long)


_closure_cache = {}

cdef bint trystrat(int n, int obs, strat, substrat, dec_nodes, np,
        int not_dec_nodes, downstream, iobs):
    """Try a specific application of Theorem 26, namely steps 2-3 in
    Appendix C

    n -- number of nodes
    obs -- bitstring specifying which nodes are observable
    strat -- choices of root notes, in same order as dec_nodes
    substrat -- ordering for strat and dec_nodes we want to attempt
    dec_nodes -- "tricky nodes", in same order as strat
    np -- bitstring of parents for each node with step 1 already applied
    not_dec_nodes -- bitstring of observable nodes not in dec_nodes
    downstream -- for each node X, bitstring of observable nodes Y with
        X squiggly arrow to Y, i.e. path from X to Y through unobserved
    iobs -- observable conditional independences for this GDAG

    In the notation of Appendix C:
        The node number of T_i is dec_nodes[substrat[i]]
        The node number of R_i is strat[substrat[i]]

    Returns True if application has succeded and so C=I for this GDAG
    """
    cdef int targ, targbit, calc, consider, x, xind, ipar, ind, i, ibit
    s_np = np[:]
    cdef int already = 0

    # Step 2
    for i in substrat:
        targ = dec_nodes[i]
        targbit = 1<<targ
        calc = strat[i]

        s_np[targ] &= not_dec_nodes | already # Step 2(a)
        for j in bitsof(downstream[calc] & ~already & ~targbit):
            if s_np[targ] & ~s_np[j] == 0:
                s_np[j] |= targbit # Step 2(b)
        already |= targbit

    # [ Step 3 is implicit - we simply ignore the unoberved nodes ]

    # Calculate the standard generating set of conditional independences
    # (node indep. non-descendants given parents for each node) and see
    # if they are all observable independences for the original GDAG
    nondesc = [obs & ~(1<<x) for x in range(n)]
    for i in bitsof(obs):
        ibit = 1<<i
        consider = s_np[i]
        while consider:
            x = consider & ~(consider - 1)
            consider &= ~x
            xind = getindex(x)
            consider |= s_np[xind]
            nondesc[xind] &= ~ibit

    for i in bitsof(obs):
        ipar = s_np[i]
        ind = nondesc[i] & ~ipar
        if ind and Triplet(1<<i, ind, ipar) not in iobs:
            return False

    return True

def isclassical(int n, par, int obs):
    """ Determine if our sufficient condition tells us C = I

    (n, par, obs) is the standard format for a GDAG:
    n -- number of nodes
    par[i] -- bitstring of parents of node i
    obs -- bitstring of observable nodes
    """

    # closure(dag_indeps(par)) is independent of obs, and many GDAGs
    # share a parent structure, so cache it by par
    key = tuple(par)
    cdag = _closure_cache.get(key)
    if cdag is None:
        if len(_closure_cache) >= 4096:
            _closure_cache.clear()
        cdag = closure(dag_indeps(key))
        _closure_cache[key] = cdag
    iobs = observable_indeps(cdag, obs)

    cdef int i, ibit, u_roots, ipar, i_is_obs, consider
    cdef int x, xind, xpar, xpar_u
    
    # Identify "tricky nodes" and the possible associated "root nodes"
    np = list(par)
    dec_nodes = []
    dec_options = []
    knows = []
    downstream = [0 for x in range(n)]
    not_dec_nodes = 0
    for i in range(n):
        ibit = 1<<i
        u_roots = 0
        ipar = par[i]
        i_is_obs = ibit & obs

        # Explore unobserved parents, and their unobserved parents etc
        consider = ipar & ~obs
        while consider:
            x = consider & ~(consider - 1)
            xind = getindex(x)
            xpar = par[xind]
            xpar_u = xpar & ~obs
            consider &= ~x
            consider |= xpar_u

            if i_is_obs:
                downstream[xind] |= ibit
            if xpar_u == 0:
                u_roots |= x
            # Step 1 of Appendix C can be applied once and for all
            np[i] |= xpar & obs

        if i_is_obs:
            if u_roots:
                dec_nodes.append(i)
                dec_options.append(bitsof(u_roots))
            else:
                not_dec_nodes |= ibit

    # If all observed nodes have only observed parents, trivially C=I
    if dec_nodes == []:
        return True

    # Now the two things we need to search over:
    # "every possible ordering of T"
    substrats = list(itertools.permutations(range(len(dec_nodes))))
    # "each element T_i associated with every possible R_i"
    for strat in itertools.product(*dec_options):
        for substrat in substrats:
            if trystrat(n, obs, strat, substrat, dec_nodes, np,
                    not_dec_nodes, downstream, iobs):
                return True

    return False

cdef inline bint implies(implicants, implicand):
    """Does (c.i Triple) implicand follow from one of the implicants?"""
    for have in implicants:
        if implicand < have:
            return True
    return False

def observable_indeps(gen, int obs):
    """Calculate all the conditional independences that follow from gen
    and only mention nodes in the bitstring obs"""
    res = []

    cdef int i, j, k
    for i in nonemptysubsetsof(obs):
        for j in nonemptysubsetsof(obs & ~i):
            for k in subsetsof(obs & ~i & ~j):
                ijk = Triplet(i,j,k)
                if implies(gen, ijk):
                    res.append(ijk)

    return tuple(res)

def dag_indeps(par):
    """Calculate a generating set of conditional independces for a DAG

    par[i] -- bit array giving parents of i, assumed to all be before i
    """
    strat = []
    for i in range(len(par)):
        C = par[i]
        if (1<<i)-1 & ~C:
            strat.append(Triplet(1<<i, (1<<i)-1 & ~C, C))
    return strat;

cdef inline int applyperm(perm, int x):
    """Apply the permutation perm to the bitstring x"""
    cdef int y = 0, i, j
    for i,j in enumerate(perm):
        if ((1<<j) & x):
            y |= (1<<i)
    return y

def irreducible(n, par, obs):
    """Is (n, par, obs) worth listing after considernig reducibility?
    
    Returns False if is (n, par, obs) reducible to a smaller GDAG for
    which either
    (a) C and I are unchanged (i.e. reduction not by 1, 4, and 5 of
        Appendix D.1), or
    (b) Our sufficient condition fails for the smaller GDAG, and so the
        smaller GDAG will already be on the list of "interesting" ones
    """

    # Appendix D.1, 4: If there is a proof that C != I without some
    # obvervable node, then making that node trivial gives a proof for
    # par
    for i in bitsof(obs):
        noti = ~(1<<i)
        thispar = [x & noti for x in par]
        thispar[i] = 0
        if not isclassical(n, thispar, obs): # case (b) above
            return False

    # Appendix D.1, 5: Causal conditionals for nodes whose parents are
    # all observable can be seen in statistics. Hence we can enforce the 
    # uselessness" of an incoming edge and then use a proof without that
    # edge
    all_p_obs = 0
    for i in bitsof(obs):
        if par[i] & ~obs == 0:
            all_p_obs |= (1<<i)

    thispar = list(par)
    for i in bitsof(all_p_obs):
        for j in propersubsetsof(par[i]):
            thispar[i] = j
            if not isclassical(n, thispar, obs): # case (b) above
                return False
        thispar[i] = par[i]

    # Appendix D.2, 1: "Heisenberg picture" / composition of channels:
    # if a node has an unobserved parent with no other children, then
    # the channel can be incorporated into the measurement / composed
    # with channel
    for i in range(n):
        for j in bitsof(par[i] & ~obs):
            jbit = 1<<j
            if any([(jbit & parents and node != i) for 
                    node, parents in enumerate(par)]):
                continue
            return False # case (a) above

    # Appendix D.2, 2: If we have a parentless unobserved node with only
    # two children, one of which is observed and has no other parents
    # then all that unobserved node can achieve is deciding that child
    # might just use that child directly instead, giving smaller DAG
    u_nodes = ((1<<n) - 1) & ~obs
    for i in bitsof(u_nodes):
        if par[i]:
            continue
        ibit = 1<<i
        children = [j for j in range(n) if (par[j] & ibit)]
        if len(children) != 2:
            continue
        for j in children:
            if (1<<j) & obs and par[j] == ibit:
                return False # case (a) above

    # Appendix D.1, 6: Unobserved nodes whos parents and children are
    # subsets of another can be subsumed into that
    for i in bitsof(u_nodes):
        ibit = 1<<i
        for j in bitsof(u_nodes):
            if i == j:
                continue

            if par[i] & ~par[j]:
                continue

            jbit = 1<<j
            kidsok = True
            for k in range(n):
                if k != j and par[k] & ibit and not par[k] & jbit:
                    kidsok = False
                    continue

            if kidsok:
                return False # case (a) above

    return True

def trydag(args):
    """Return args=(n, par, obs) if an irreducibily interesting GDAG"""
    if not isclassical(*args) and irreducible(*args):
        return args

def trydag_reducible(args):
    """Does args=(n, par, obs) give an interesting GDAG?"""
    if not isclassical(*args):
        return args

def trydag_count(args):
    """Does args=(n, par, obs) give an interesting GDAG?"""
    return isclassical(*args)

def find_candidates(int n, bint irrcheck):
    """Returns a list of GDAGs with n nodes that are worth looking at

    This means that:
    (a) No GDAG in the list is isomorphic to another
    (b) If irrcheck is True, all GDAGs pass initial reducibilty checks:
        (i)   They are connected
        (ii)  They don't have childless unobserved nodes
        (iii) They don't have unobserved nodes with a sole parent,
              also unobserved

    It is worth checking for a few reducibilities now since it avoids
    expensive isomorphism tests for obviously uninteresting GDAGs

    A GDAG is packed into a single integer key: bit n + j*(j-1)/2 + k
    set means k is a parent of j (the layout the loop over i already
    enumerates), and the low n bits are obs. done is a flat bitmap
    over keys, so the duplicate check is one bit test. When a new
    GDAG is found, every permuted form of it is marked in done, using
    precomputed tables of where each key bit moves under each
    permutation (or -1 if the permuted edge would run from a higher
    to a lower node, which no enumerated key can match anyway)"""
    cdef int numbits = n*(n-1) // 2
    cdef int everything = (1<<n) - 1
    cdef int totalbits = numbits + n
    cdef unsigned long long maxkey = 1ULL << totalbits
    cdef int maxi = 1<<numbits
    cdef int i, j, c, p, b, obs, tgt, pi, bit, ok
    cdef int childless = 0, connectedA, oldconnectedA, ipar, inodes
    cdef unsigned long long key, src, out
    cdef int sigma[32]
    cdef int par_c[32]

    perms = list(itertools.permutations(range(n)))
    cdef int nperms = len(perms)
    cdef int *perm_map = <int *> malloc(nperms * totalbits * sizeof(int))
    cdef unsigned char *done = \
        <unsigned char *> malloc(<size_t> (maxkey // 8 + 1))
    if perm_map == NULL or done == NULL:
        free(perm_map)
        free(done)
        raise MemoryError()
    memset(done, 0, <size_t> (maxkey // 8 + 1))

    for pi, perm in enumerate(perms):
        # perm[new] = old, so sigma maps old node numbers to new
        for j in range(n):
            sigma[perm[j]] = j
        for j in range(n): # observability bits
            perm_map[pi*totalbits + j] = sigma[j]
        for c in range(n): # edge bits: parent p of child c
            for p in range(c):
                b = n + c*(c-1)//2 + p
                if sigma[p] < sigma[c]:
                    perm_map[pi*totalbits + b] = \
                        n + sigma[c]*(sigma[c]-1)//2 + sigma[p]
                else:
                    perm_map[pi*totalbits + b] = -1

    trythese = set()

    try:
        for i in range(maxi):
            for j in range(n):
                par_c[j] = (i >> (j*(j-1) // 2)) & ((1<<j) - 1)

            if irrcheck:
                childless = everything
                for j in range(n):
                    childless &= ~par_c[j]

                # Appendix D.1, 1: disconnected DAGs don't tell us more
                # than each component does
                connectedA = 1
                oldconnectedA = 0
                while connectedA != oldconnectedA:
                    oldconnectedA = connectedA
                    for j in range(n):
                        inodes = par_c[j] | (1<<j)
                        if inodes & connectedA:
                            connectedA |= inodes
                if connectedA != everything:
                    continue

            par = None
            for obs in range(1<<n):
                if irrcheck:
                    # Appendix D.1, 2: Childless unobserved nodes are
                    # pointless
                    if childless & ~obs:
                        continue

                    # Appendix D.1, 3: "Schrodinger picture": if an
                    # unobserved node only has one parent, and its
                    # unobserved, just apply channel to parent
                    ok = 1
                    for j in range(n):
                        if not ((1<<j) & obs):
                            ipar = par_c[j]
                            if ipar & (ipar - 1) == 0 and ipar & ~obs:
                                ok = 0
                                break
                    if not ok:
                        continue

                key = (<unsigned long long> i << n) | <unsigned int> obs
                if done[key >> 3] & (1 << (key & 7)):
                    continue

                # New isomorphism class: mark every permuted form
                for pi in range(nperms):
                    src = key
                    out = 0
                    ok = 1
                    while src:
                        bit = _ctz(src)
                        src &= src - 1
                        tgt = perm_map[pi*totalbits + bit]
                        if tgt < 0:
                            ok = 0
                            break
                        out |= 1ULL << tgt
                    if ok:
                        done[out >> 3] |= 1 << (out & 7)

                if par is None:
                    par = tuple([par_c[j] for j in range(n)])
                trythese.add((n, par, obs))
    finally:
        free(perm_map)
        free(done)

    return trythese

def enumdags(int n, int numprocesses):
    """Main interface to this module: output a list of interesting
    GDAGs of size n, using numprocesses processes in parallel"""
    trythese = find_candidates(n, True)
    pool = multiprocessing.Pool(processes=numprocesses)
    # sorted + chunksize keeps same-par GDAGs in the same worker, so
    # each worker's closure cache gets hits
    res = pool.imap_unordered(trydag, sorted(trythese), chunksize=64)
    for x in res:
        if x is not None:
            (n, par, obs) = x
            print([val2nice(x) for x in par], val2nice(obs))

def enumdags_reducible(int n, int numprocesses):
    """Main interface to this module: output a list of interesting
    GDAGs of size n, using numprocesses processes in parallel"""
    trythese = find_candidates(n, False)
    pool = multiprocessing.Pool(processes=numprocesses)
    res = pool.imap_unordered(trydag_reducible, sorted(trythese),
                              chunksize=64)
    for x in res:
        if x is not None:
            (n, par, obs) = x
            print([val2nice(x) for x in par], val2nice(obs))

def countdags(int n, int numprocesses):
    """Alternative interface to this module: output n, the number of
    GDAGs, and the number of interesting ones, using numprocesses
    processes in parallel"""
    trythese = find_candidates(n, False)
    pool = multiprocessing.Pool(processes=numprocesses)
    res = pool.imap_unordered(trydag_count, sorted(trythese),
                              chunksize=64)
    print(n, len(trythese), sum(res))
