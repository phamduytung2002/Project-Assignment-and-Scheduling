from IP.IP import ILP
from CP.CP import CP
from option.options import args_parser
if __name__ == '__main__':
    args = args_parser()
    path = 'data/'+args.dataset
    if args.method == 'IP':
        ILP(path)
    if args.method == 'CP':
        CP(path)
