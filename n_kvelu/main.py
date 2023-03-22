from ctool import OpCount


def main():
    OpCount.op("add", "3")
    OpCount.op("add", "3")
    OpCount.op("add", "3")
    OpCount.op("add", "3")
    OpCount.print_results()
    OpCount.clean()
    OpCount.print_results()


if __name__ == "__main__":
    main()
