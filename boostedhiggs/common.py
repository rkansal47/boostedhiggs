import numpy as np
def getParticles(events,id=25,flags=['fromHardProcess', 'isLastCopy'],status=None):
    absid = np.abs(events.GenPart.pdgId)
    idx = (absid == id) & events.GenPart.hasFlags(['fromHardProcess', 'isLastCopy'])
    if status:
        idx = idx & events.GenPart.status == status
    return events.GenPart[idx],idx
