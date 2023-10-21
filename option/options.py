import argparse
def args_parser():
    parser = argparse.ArgumentParser()
    # save file 
    parser.add_argument('--method', type=str, default='IP',
                        help="")
    parser.add_argument('--dataset', type=str, default='input1.txt',
                        help="data input")
    args = parser.parse_args()
    return args