import numpy as np
def getParticles(events,id=25,flags=['fromHardProcess', 'isLastCopy'],status=None):
    absid = np.abs(events.GenPart.pdgId)
    if isinstance(id,int):
        idx = (absid == id) & events.GenPart.hasFlags(flags)
    else:
        idx = (absid >= id[0]) & (absid <= id[1]) & events.GenPart.hasFlags(flags)
    if status:
        idx = idx & events.GenPart.status == status
    return events.GenPart[idx],idx

def match(left, right, metric, maximum=np.inf):
    '''Matching utility                                                                                                                                                         
    For each item in ``left``, find closest item in ``right``, using function ``metric``.                                                                                      
    The function must accept two broadcast-compatible arrays and return a numeric array.                                                                                         
    If maximum is specified, mask matched elements where metric was greater than it.                                                                                               
    '''
    lr = left.cross(right, nested=True)
    mval = metric(lr.i0, lr.i1)
    idx = mval.argmin()
    print(mval,idx)
    if maximum < np.inf:
        matched = lr.i1[idx[mval[idx] < maximum]]
        return matched.copy(content=matched.content.pad(1)).flatten(axis=1)
    else:
        return lr.i1[idx]
    
def getFlavor(childid, momid=None):
    x = ((childid == 24).any() * 6 + (childid == 13).any() * 5 + (childid == 11).any() * 4 + (childid == 5).any() * 3 + (childid == 4).any() * 2 + (childid < 4).all() * 1).pad(1, clip=True).fillna(0).flatten()
    x = x + ((momid == 5).any() * 7)
    return x
    
def matchedParticleFlavor(candidates, particles, types='child', maxdR=0.8): 
    matched = match(candidates, particles, lambda a, b: a.delta_r(b), maxdR)
    if types=='child':
        retid = abs(matched.children.pdgId)
    else:
        retid = abs(matched.pdgId)
    return retid

