import os
import sys
import glob


def which(program):
    """ returns the path to an executable or None if it can't be found"""
    def is_exe(_fpath):
        return os.path.isfile(_fpath) and os.access(_fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ['PATH'].split(os.pathsep):
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file
    return None


def find_markers_file(opts, ext='.bed'):
    libs_dir = os.path.dirname(__file__)
    markers_dir = os.path.join(libs_dir, 'markers')

    if opts.markers:
        markers_file = opts.markers
    else:
        genome = opts.genome
        if genome in ['GRCh37', 'hg19']:
            markers_file = os.path.join(markers_dir, 'GRCh37.autosomes.phase3_shapeit2_mvncall_integrated.20130502.SNV.genotype.sselect_v4_MAF_0.4_LD_0.8' + ext)
        elif genome in ['GRCh38', 'hg38']:
            markers_file = os.path.join(markers_dir, 'GRCh38.autosomes.phase3_shapeit2_mvncall_integrated.20130502.SNV.genotype.sselect_v4_MAF_0.4_LD_0.8' + ext)
        elif genome in ['GRCm38', 'mm10', 'mouse']:
            markers_file = os.path.join(markers_dir, 'SureSelect_Mouse_All_Exon_V1_GRCm38_markers_MAF_0.4_LD_0.8' + ext)
        else:
            print('Genome ' + opts.genome + ' is not recognized. Available: GRCh37, GRCh38, GRCm38')
            sys.exit(2)

    if not os.path.exists(markers_file):
        print('ERROR: Marker file {0} cannot be found.'.format(markers_file))
        sys.exit(2)

    return markers_file

