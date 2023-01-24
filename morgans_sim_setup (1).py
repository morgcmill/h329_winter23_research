import pynbody


def sim_setup(path_to_sim):
    """It's always good to comment your codes for later!"""
    # load in the simulation
    print("loading sim")
    sim = pynbody.load(path_to_sim)
    # put it in the right units
    print("doing units")
    sim.physical_units()
    # select the main halo?
    print("selecting halos")
    h = sim.halos()
    h1 = h[1]
    # center it or turn it face on
    print("making faceon")
    pynbody.analysis.angmom.faceon(h1)
    
    return sim, h, h1