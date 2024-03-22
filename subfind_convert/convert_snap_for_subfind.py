#!/bin/env python

import os
import h5py


def main():
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('--snap', type=int)
    parser.add_argument('--directory', type=str)
    parser.add_argument('-v', '--verbose', action='count', default=0)
    args = parser.parse_args()

    verbose = args.verbose > 0

    snap_number = int(args.snap)
    directory = args.directory

    print(verbose, snap_number, directory)

    convert(snap_number, directory, verbose)

def convert(snap_number, in_directory, out_directory=None, verbose=False):
    if out_directory is None:
        out_directory = in_directory

    snap = f'snapshot_{snap_number:04d}.hdf5'
    out_snap = f'subsnap_{snap_number:03d}.hdf5'

    if verbose:
        print((f'Input Directory: {in_directory}\n Output Directory: {out_directory}'
               f'Snap number: {snap_number}\nOutput snap name: {out_snap}\n'))

    output_units = {
        'UnitLength_in_cm': 3.085678e+21,
        'UnitMass_in_g': 1.989e+43,
        'UnitVelocity_in_cm_per_s': 100000.0, 
    }

    attrs_update = [
        'CanHaveTypes',
        'Flag_Entropy_ICs',
        'InitialMassTable',
        'MassTable',
        'NumPart_ThisFile',
        'NumPart_Total',
        'NumPart_Total_HighWord',
        'TotalNumberOfParticles',
    ]

    out_filename = f'{out_directory}/{out_snap}'
    if os.path.exists(out_filename):
        if verbose:
            print("Deleting previous output file...")

        os.remove(out_filename)

    with h5py.File(out_filename, 'a') as hf_out:

        hf_out.create_group('Header')

        with h5py.File(f'{in_directory}/{snap}', 'r') as hf_in:

            scale_factor = hf_in['Cosmology'].attrs['Scale-factor']
            hubble_param = hf_in['Cosmology'].attrs['H [internal units]'] / 100

            fields = [
                ['PartType0/ParticleIDs', 'PartType0/ParticleIDs', None],
                ['PartType0/Coordinates', 'PartType0/Coordinates',
                    output_units['UnitLength_in_cm'] * scale_factor / hubble_param],
                ['PartType0/Masses', 'PartType0/Masses',
                    output_units['UnitMass_in_g'] / hubble_param],
                ['PartType0/InternalEnergies', 'PartType0/InternalEnergy',
                    output_units['UnitVelocity_in_cm_per_s']**2],
                ['PartType0/Densities', 'PartType0/Density',
                    output_units['UnitMass_in_g'] * output_units['UnitLength_in_cm']**-3 *
                            scale_factor**-3 * hubble_param**2],
                ['PartType0/Velocities', 'PartType0/Velocities',
                    output_units['UnitVelocity_in_cm_per_s'] * scale_factor**0.5],
                ['PartType0/SmoothingLengths', 'PartType0/SubfindHsml',
                    output_units['UnitLength_in_cm'] * scale_factor / hubble_param],

                ['PartType1/ParticleIDs', 'PartType1/ParticleIDs', None],
                ['PartType1/Coordinates', 'PartType1/Coordinates',
                    output_units['UnitLength_in_cm'] * scale_factor / hubble_param],
                ['PartType1/Masses', 'PartType1/Masses',
                    output_units['UnitMass_in_g'] / hubble_param],
                ['PartType1/Velocities', 'PartType1/Velocities',
                    output_units['UnitVelocity_in_cm_per_s'] * scale_factor**0.5],

                ['PartType4/ParticleIDs', 'PartType4/ParticleIDs', None],
                ['PartType4/Coordinates', 'PartType4/Coordinates',
                    output_units['UnitLength_in_cm'] * scale_factor / hubble_param],
                ['PartType4/Masses', 'PartType4/Masses',
                    output_units['UnitMass_in_g'] / hubble_param],
                ['PartType4/Velocities', 'PartType4/Velocities',
                    output_units['UnitVelocity_in_cm_per_s'] * scale_factor**0.5],
                ['PartType4/SmoothingLengths', 'PartType4/SubfindHsml',
                    output_units['UnitLength_in_cm'] * scale_factor / hubble_param],

                ['PartType5/ParticleIDs', 'PartType5/ParticleIDs', None],
                ['PartType5/Coordinates', 'PartType5/Coordinates',
                    output_units['UnitLength_in_cm'] * scale_factor / hubble_param],
                ['PartType5/SubgridMasses', 'PartType5/Masses',
                    output_units['UnitMass_in_g'] / hubble_param],
                ['PartType5/Velocities', 'PartType5/Velocities',
                    output_units['UnitVelocity_in_cm_per_s'] * scale_factor**0.5],
                ['PartType5/SmoothingLengths', 'PartType5/SubfindHsml',
                    output_units['UnitLength_in_cm'] * scale_factor / hubble_param],
            ]

            ## set up input units
            input_units = {
                'UnitLength_in_cm': 3.08568e+24,
                'UnitMass_in_g': 1.98841e+43,
                'UnitTime_in_s': 3.08568e+19,

            }
            input_units['UnitVelocity_in_cm_per_s'] = input_units['UnitLength_in_cm'] /\
                    input_units['UnitTime_in_s']

            hf_out['Header'].attrs.update(hf_in['Header'].attrs)

            ## add cosmo information
            hf_out['Header'].attrs['Omega0'] = 0.3
            hf_out['Header'].attrs['OmegaLambda'] = 0.7
            hf_out['Header'].attrs['HubbleParam'] = 0.6711

            ## add extra flags
            hf_out['Header'].attrs['Flag_Cooling'] = 0
            hf_out['Header'].attrs['Flag_DoublePrecision'] = 0
            hf_out['Header'].attrs['Flag_Feedback'] = 0
            hf_out['Header'].attrs['Flag_IC_Info'] = 0
            hf_out['Header'].attrs['Flag_Metals'] = 11
            hf_out['Header'].attrs['Flag_Sfr'] = 0
            hf_out['Header'].attrs['Flag_StellarAge'] = 0

            ## fix boxsize dimensions (3 -> 1)
            hf_out['Header'].attrs['BoxSize'] = hf_out['Header'].attrs['BoxSize'][0]

            ## Remove extra PartType (7 -> 6)
            hf_out['Header'].attrs['NumPartTypes'] = 6

            for attr in attrs_update:
                hf_out['Header'].attrs[attr] = hf_out['Header'].attrs[attr][:-1]


            for field in fields:
                if verbose:
                    print(f'Writing {field}')

                temp = hf_in[field[0]][:]

                ## get SWIFT conversions to CGS
                conv_factor = hf_in[field[0]].attrs[('Conversion factor to physical CGS '
                                                     '(including cosmological corrections)')] 

                if conv_factor != 1.:

                    conv_factor /= field[2]
                    temp *= conv_factor

                    ## Convert CGS to GADGET units
                    hf_out[field[1]] = temp #/ field[2]
                else:
                    hf_out[field[1]] = temp

if __name__ == "__main__":
    main()
