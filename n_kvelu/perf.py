


class PerfClass(object):

    def __new__(cls):
    if not hasattr(cls, 'instance'):
      cls.instance = super(PerfClass, cls).__new__(cls)
    cls.nr_add = 0 # Number of Additions
    cls.nr_mult = 0 # Number of Multiplications
    cls.nr_inv = 0 # Number of Inversion
    cls.nr_sq = 0 # Number of Squarings
    cls.nr_frob = 0 # Number of frobenius_endomorphism

    return cls.instance


    def print_counters():
        counter_perf = PerfClass()

        print("-- \t Results \t--")
        print("Additions:\t {}".format(counter_perf.nr_add))
