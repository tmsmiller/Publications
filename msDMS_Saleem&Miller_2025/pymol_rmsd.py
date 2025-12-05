from pymol import cmd
import argparse
from glob import glob as gb

def get_models(prefix, dirpath):
    models = []
    for model in gb(dirpath + prefix + '*'):
        if model.split('.')[-1] != 'cif':
            continue
        models.append(model)
    #print(models)
    return models

def get_rmsd(reference_obj, predicted):
    pred_obj = predicted.split('\\')[-1].split('.')[0]

    print(reference_obj, pred_obj)
    cmd.load(predicted, pred_obj)
    #cmd.select(reference_obj)
    #cmd.select(pred_obj)
    out = cmd.align(reference_obj, pred_obj, cycles=0)
    cmd.delete(pred_obj)
    print(out[0], out[3])
    return out[3]

def write_rmsds(rmsd_dict, outname):
    out = open(outname, 'w')
    for model in rmsd_dict.keys():
        out.write(model.split('\\')[-1].split('.')[0] + '\t' + str(rmsd_dict[model]) + '\n')

def _rmsds_main(ref_struct, pred_prefix, output, directory):
    model_list = get_models(pred_prefix, directory)
    ref_obj = ref_struct.split('/')[-1].split('.')[0]
    cmd.load(ref_struct, ref_obj)
    rmsds = {}
    for model in model_list:
        rmsds[model] = get_rmsd(ref_obj, model)
    write_rmsds(rmsds, output)

if __name__ == '__main__':
    argparser = argparse.ArgumentParser(description='Compare many structure predictions to a single pdb structure. Output is RMSD in angstroms for all atoms')
    argparser.add_argument('--pred_dir', default='./', help='directory to search for predicted structures')
    argparser.add_argument('--pred_prefix', required=True, help='prefix for all predicted structures to compare')
    argparser.add_argument('--ref_struct', required=True, help='reference structure to compare against')
    argparser.add_argument('--output', required=True, help='output filename including file ext e.g. (output.txt)')
    args = argparser.parse_args()

    _rmsds_main(args.ref_struct, args.pred_prefix, args.output, args.pred_dir)
    



