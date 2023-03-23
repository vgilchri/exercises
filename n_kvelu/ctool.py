import matplotlib.pyplot as plt
class OpCount(object):
    nr_add = 0
    nr_mult = 0
    nr_sq =  0
    nr_sqr = 0
    nr_inv = 0

    field_op = {}

    __instance = None

    def __new__(cls, *args, **kwargs):
        if not OpCount.__instance:
            OpCount.__instance = object.__new__(cls)
        return OpCount.__instance

    @staticmethod
    def clean():
        del OpCount.field_op
        OpCount.field_op = {}

    @staticmethod
    def op(operation, field=None):
        if field == None:
            field = "1"
        try:
            OpCount.field_op[field][operation] += 1
        except:
            try:
                t = OpCount.field_op[field]
                try:
                    a = OpCount.field_op[field][operation]
                except:
                    OpCount.field_op[field][operation] = 1
            except:
                OpCount.field_op[field] = {}
                OpCount.field_op[field][operation] = 1

    @staticmethod
    def print_results():
        print("-- \t Printing number of Operations \t --")
        for k in OpCount.field_op:
            print("-- \t Operations in k = {} \t --".format(k))
            for op in OpCount.field_op[k]:
                print("\t op: {} - nr: {}".format(op, OpCount.field_op[k][op]))
        print("---------------------------------------------")


    @staticmethod
    def plot(file_name):
        for k in OpCount.field_op:
            names = list(OpCount.field_op[k].keys())
            values = list(OpCount.field_op[k].values())
            plt.bar(range(len(OpCount.field_op[k])), values, tick_label=names)
            file_name = file_name + "_" + str(k) + ".png"
            plt.savefig(file_name)
