from aflow import keywords
import aflow
import decimal
from decimal import Decimal
from pymatgen import MPRester
from tabulate import tabulate
from itertools import permutations, accumulate
import os
import numpy as np
import pymatgen.analysis.structure_matcher
from pymatgen.core import Structure
from sys import argv
import smtplib
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
from email.mime.base import MIMEBase
from email import encoders
file1= str(argv[1])
MATCOR_out = str(argv[2])
materials_in = []
with open(file1, 'r') as infile:
    for inf in infile:
        materials_in.append(inf)
materials = []
for mot in materials_in:
    materials.append(mot.replace('\n', ''))

# edit
mpr = MPRester('YOUR API ID')
base_property_MP_1 = 'density'
base_property_MP_2 = ''
base_property_AFLOW = 'density'
max_difference = 0.0
ltol_in = 0.2
stol_in = 0.3
angle_tol_in = 5



key_headers = 0
key_MATCOR = ''
supported_properties_py = mpr.supported_properties


def gcd(a, b):
    while b > 0:
        a, b = b, a % b
    return a


properties1 = []
body3 = []


def pymatgen_to_AFLOW(id, mpr, base_property_py1, base_property_py2, base_property_af, ltol, stol, angle_tol):
    # avoids possible reference errors
    auid_AFLOW_base = ''
    auid_reduced_and_hubbard = ''
    final_sg_aflow = ''
    space_group_pymatgen = ''
    all_hubbards_aflow_reduced_and_hubbard = ''
    formula_input = ''
    is_structure_match = ''
    structrure_verification = ''
    species_CONT1 = ''
    composition_CONT = ''
    is_hubbard_py = ''
    hubbard_U = ''
    error_py_AFLOW = ''
    # Properties for pymatgen
    criteria1 = {'material_id': id}
    properties1 = supported_properties_py
    result_pymatgen = mpr.query(criteria=criteria1, properties=properties1)
    if len(result_pymatgen) != 0:
        input_structure = mpr.get_structure_by_material_id(id)
        # defines base property for MP
        for elem in result_pymatgen:
            property_py = elem.get(base_property_py1)
            is_hubbard_py = elem.get('is_hubbard')
            formula_input = str(elem.get('pretty_formula'))
            if len(base_property_py2) == 0:
                if property_py is None:
                    base_property_py1 = 'None'
                else:
                    base_property_py1 = elem.get(base_property_py1)
                break
            if len(base_property_py2) != 0:
                if property_py is None:
                    base_property_py1 = 'None'
                else:
                    base_property_py1 = elem.get(base_property_py1).get(base_property_py2)
        # gets pretty formula
        for prop in result_pymatgen:
            substance_py = prop.get('pretty_formula')
            elements1 = (prop.get('elements'))
            elements1.sort()
            space_group = prop.get('spacegroup')
            space_group_pymatgen = str(space_group.get('symbol')) + ' #' + str(space_group.get('number'))

            length_loop = []
            all_hubbards_aflow = []
            all_auid = []
            all_compostion = []
            all_loops = []
            all_space_groups = []
            all_base_property = []
            Kn = 1
            length_elements1 = len(elements1)
            K = (keywords.species % elements1[0])
            while (length_elements1 > 1) & (Kn < length_elements1):
                K = K & (keywords.species % elements1[Kn])
                Kn = Kn + 1

            K = K & (keywords.nspecies == length_elements1)

            # add specific catalog
            result_aflow_rep = aflow.search().filter(K)

            all_space_groups_0 = []
            try:
                for entryrep in result_aflow_rep:
                    all_base_property.append(getattr(entryrep, base_property_af))
                    all_hubbards_aflow.append(entryrep.ldau_TLUJ)
                    all_compostion.append(entryrep.composition)
                    all_auid.append(entryrep.auid)
                    all_loops.append(entryrep.loop)
                    all_space_groups.append(entryrep.sg)
                    all_space_groups_0.append(entryrep.sg2)
            except AssertionError:
                auid_AFLOW_base = 'NO_MATCH'
                error_py_AFLOW = 'NO RESPONSE FROM AFLOW'

            def sg_representative(all_space_groups1):
                # filters through space group
                all_hubbards_aflow_reduced_and_sg = []
                auid_reduced_and_sg = []
                base_property_and_sg = []
                composition_reduced_and_sg = []
                loop_reduced_and_sg = []
                all_reduced_and_sg_pos = []

                all_space_groups_and_sg = []
                edited_space_groups_AFLOW = []
                length_sg = []
                edited_sg0_2 = []
                length_sg0_2 = []
                for space_groups in all_space_groups1:
                    length_sg.append(len(space_groups))
                    for space_group_AFLOW in space_groups:
                        space_group_AFLOW1 = space_group_AFLOW.replace('{', '')
                        space_group_AFLOW1 = space_group_AFLOW1.replace('}', '')
                        space_group_AFLOW1 = space_group_AFLOW1.replace('\n', '')
                        edited_space_groups_AFLOW.append(space_group_AFLOW1)
                edited_all_space_groups = [edited_space_groups_AFLOW[x - y: x] for x, y in zip(
                    accumulate(length_sg), length_sg)]

                for edit in edited_all_space_groups:
                    length_sg0_2.append(len(edit) - 1)
                    edited_sg0_2.append(edit[0])
                    if len(edit) >= 3:
                        edited_sg0_2.append(edit[2])
                    else:
                        edited_sg0_2.append(edit[0])
                edited_all_sg0_2 = [edited_sg0_2[x - y: x] for x, y in zip(
                    accumulate(length_sg0_2), length_sg0_2)]

                for space_group_AFLOW in edited_all_sg0_2:
                    if space_group_pymatgen in space_group_AFLOW:
                        # makes list of all who do have same space group, formula
                        all_reduced_and_sg_pos = [i for i, j in enumerate(edited_all_sg0_2) if
                                                  space_group_pymatgen in j]
                    if space_group_pymatgen not in space_group_AFLOW:
                        continue

                if len(all_reduced_and_sg_pos) != 0:
                    for reduced_and_sg_pos in all_reduced_and_sg_pos:
                        all_hubbards_aflow_reduced_and_sg.append(all_hubbards_aflow[reduced_and_sg_pos])
                        auid_reduced_and_sg.append(all_auid[reduced_and_sg_pos])
                        composition_reduced_and_sg.append(all_compostion[reduced_and_sg_pos])
                        base_property_and_sg.append(all_base_property[reduced_and_sg_pos])
                        loop_reduced_and_sg.append(all_loops[reduced_and_sg_pos])
                        all_space_groups_and_sg.append(all_space_groups1[reduced_and_sg_pos])
                if len(all_reduced_and_sg_pos) == 0:
                    auid_reduced_and_sg = 'place_holder'
                return [all_hubbards_aflow_reduced_and_sg, auid_reduced_and_sg, composition_reduced_and_sg,
                        base_property_and_sg, loop_reduced_and_sg, all_space_groups_and_sg]

            x = sg_representative(all_space_groups1=all_space_groups)
            all_hubbards_aflow_reduced_and_sg = x[0]
            auid_reduced_and_sg = x[1]
            composition_reduced_and_sg = x[2]
            base_property_and_sg = x[3]
            loop_reduced_and_sg = x[4]
            all_space_groups_and_sg = x[5]

            if auid_reduced_and_sg == 'place_holder':
                y = sg_representative(all_space_groups1=all_space_groups_0)
                all_hubbards_aflow_reduced_and_sg = y[0]
                auid_reduced_and_sg = y[1]
                composition_reduced_and_sg = y[2]
                base_property_and_sg = y[3]
                loop_reduced_and_sg = y[4]
                all_space_groups_and_sg = y[5]

            # filters through pretty formula
            all_hubbards_aflow_reduced_and_subs = []
            auid_reduced_and_subs = []
            composition_reduced_and_subs = []
            base_property_and_subs = []
            loop_reduced_and_subs = []
            all_edited_substances_AFLOW = []
            all_space_groups_and_subs = []
            for compo in composition_reduced_and_sg:
                atoms = compo
                # changes AFLOW to pymatgen format
                result = atoms[0]
                for i in atoms[1:]:
                    result = gcd(result, i)
                pretty_natoms = []
                for x in atoms:
                    pretty_natoms.append(int(x / result))
                for n, s in enumerate(pretty_natoms):
                    if s == 1:
                        pretty_natoms[n] = ' dummy'
                elements1[-1] = elements1[-1].replace('\n', '')
                substance_vector = []
                substance_vector1 = []
                for e in elements1:
                    position_atom = elements1.index(e)
                    substance_vector.append(e + str(atoms[position_atom]))
                    substance_vector1.append(e + str(pretty_natoms[position_atom]))
                substance_full = ''.join(substance_vector1)
                if 'dummy' in substance_full:
                    substance_full = substance_full.replace(' dummy', '1')
                substances_AFLOW = []
                perm = permutations(substance_vector1)
                for i in perm:
                    j = str(i)
                    j = j.replace(' dummy', '')
                    j = j.replace('(', '')
                    j = j.replace("'", '')
                    j = j.replace(')', '')
                    j = j.replace(', ', '')
                    substances_AFLOW.append(j)
                all_edited_substances_AFLOW.append(substances_AFLOW)

            all_reduced_and_substance_pos = [i for i, j in enumerate(all_edited_substances_AFLOW) if
                                             substance_py in j]
            check = []

            if len(all_reduced_and_substance_pos) != 0:
                for reduced_and_substance_pos in all_reduced_and_substance_pos:
                    check.append(all_edited_substances_AFLOW[reduced_and_substance_pos])
                    all_hubbards_aflow_reduced_and_subs.append(
                        all_hubbards_aflow_reduced_and_sg[reduced_and_substance_pos])
                    composition_reduced_and_subs.append(composition_reduced_and_sg[reduced_and_substance_pos])
                    auid_reduced_and_subs.append(auid_reduced_and_sg[reduced_and_substance_pos])
                    base_property_and_subs.append(base_property_and_sg[reduced_and_substance_pos])
                    loop_reduced_and_subs.append(loop_reduced_and_sg[reduced_and_substance_pos])
                    all_space_groups_and_subs.append(all_space_groups_and_sg[reduced_and_substance_pos])
            for c in check:
                assert substance_py in c

            # Hubbard-U filtering

            all_hubbards_aflow_reduced_and_hubbard = []
            auid_reduced_and_hubbard = []
            icsd_reduced_and_hubbard = []
            base_property_reduced_and_hubbard = []
            loop_reduced_and_hubbard = []
            all_space_groups_and_hubbard = []

            # first, if MP does not use Hubbard-U, tries to obtain all entries that don't use Hubbard U in AFLOW
            if is_hubbard_py is False:
                all_reduced_and_hubbard_pos_false = [i for i, j in enumerate(all_hubbards_aflow_reduced_and_subs) if
                                                     j is None]
                if len(all_reduced_and_hubbard_pos_false) != 0:
                    for hubbard_pos in all_reduced_and_hubbard_pos_false:
                        all_hubbards_aflow_reduced_and_hubbard.append(all_hubbards_aflow_reduced_and_subs[hubbard_pos])
                        auid_reduced_and_hubbard.append(auid_reduced_and_subs[hubbard_pos])
                        icsd_reduced_and_hubbard.append(composition_reduced_and_subs[hubbard_pos])
                        loop_reduced_and_hubbard.append(loop_reduced_and_subs[hubbard_pos])
                        base_property_reduced_and_hubbard.append(base_property_and_subs[hubbard_pos])
                        all_space_groups_and_hubbard.append(all_space_groups_and_subs[hubbard_pos])
                # if none are false, it means that the hubbard U filtering does not matter because there are no available
                # materials for which no hubbard U was used
                if len(all_reduced_and_hubbard_pos_false) == 0:
                    all_hubbards_aflow_reduced_and_hubbard = all_hubbards_aflow_reduced_and_subs
                    auid_reduced_and_hubbard = auid_reduced_and_subs
                    icsd_reduced_and_hubbard = composition_reduced_and_subs
                    loop_reduced_and_hubbard = loop_reduced_and_subs
                    base_property_reduced_and_hubbard = base_property_and_subs
                    all_space_groups_and_hubbard = all_space_groups_and_subs
            if is_hubbard_py is True:
                all_reduced_and_hubbard_pos_true = [i for i, j in enumerate(all_hubbards_aflow_reduced_and_subs) if
                                                    j is not None]

                if len(all_reduced_and_hubbard_pos_true) != 0:
                    # hubbard U used
                    for hubbard_neg in all_reduced_and_hubbard_pos_true:
                        all_hubbards_aflow_reduced_and_hubbard.append(all_hubbards_aflow_reduced_and_subs[hubbard_neg])
                        auid_reduced_and_hubbard.append(auid_reduced_and_subs[hubbard_neg])
                        icsd_reduced_and_hubbard.append(composition_reduced_and_subs[hubbard_neg])
                        loop_reduced_and_hubbard.append(loop_reduced_and_subs[hubbard_neg])
                        base_property_reduced_and_hubbard.append(base_property_and_subs[hubbard_neg])
                        all_space_groups_and_hubbard.append(all_space_groups_and_subs[hubbard_neg])
                if len(all_reduced_and_hubbard_pos_true) == 0:
                    # hubbard U not used
                    all_hubbards_aflow_reduced_and_hubbard = all_hubbards_aflow_reduced_and_subs
                    auid_reduced_and_hubbard = auid_reduced_and_subs
                    icsd_reduced_and_hubbard = composition_reduced_and_subs
                    loop_reduced_and_hubbard = loop_reduced_and_subs
                    base_property_reduced_and_hubbard = base_property_and_subs
                    all_space_groups_and_hubbard = all_space_groups_and_subs
            assert len(auid_reduced_and_hubbard) == len(icsd_reduced_and_hubbard) == len(
                loop_reduced_and_hubbard) == len(
                base_property_reduced_and_hubbard) == len(all_space_groups_and_hubbard)

            # for G_VRH error: prioritizes matching of integers, no Hubbard U filtering used
            check_base_result = [s for s in base_property_reduced_and_hubbard if s is not None]
            if len(check_base_result) == 0:
                all_hubbards_aflow_reduced_and_hubbard = all_hubbards_aflow_reduced_and_subs
                auid_reduced_and_hubbard = auid_reduced_and_subs
                base_property_reduced_and_hubbard = base_property_and_subs
                all_space_groups_and_hubbard = all_space_groups_and_subs
            print('Entries with same sg and formula, and potential hubbard U selection')
            print(auid_reduced_and_hubbard)

            # gets only initial and final space grouop for log file
            for sg_hubbard in all_space_groups_and_hubbard:
                sg_hubbard.pop(1)
            final_sg_aflow = []
            for sg2 in all_space_groups_and_hubbard:
                sg2[1] = sg2[1].replace('\n', '')
                sg = [str(element) for element in sg2]
                sg = ",".join(sg)
                final_sg_aflow.append(sg)
            min_difference_rep = []
            min_integers_rep = []
            auid_AFLOW_base = ''

            # filters through structure matching
            def structurefiltering(species_CONT, comp_cont2):
                a0 = result_aflow_CONTCAR[0]
                files = a0.files
                if 'CONTCAR.relax' in files:
                    pos_file = files.index('CONTCAR.relax')
                    MATCOR_tar = files[pos_file]('MATCOR_AFLOW_CONTCAR')
                else:
                    print('CONTCAR.relax file not provided by AFLOW')
                # makes changes for CONTCAR file in MP format
                with open('MATCOR_AFLOW_CONTCAR', 'r') as f:
                    text = f.readlines()
                    text[5] = species_CONT + '\n' + comp_cont2 + '\n'
                    if 'Direct' not in str(text[6]):
                        text.pop(7)
                my_file = open('MATCOR_AFLOW_CONTCAR', "w")
                new_file_contents = "".join(text)
                my_file.write(new_file_contents)
                my_file.close()
                output_struct = Structure.from_file('MATCOR_AFLOW_CONTCAR')
                structure_matcher = pymatgen.analysis.structure_matcher.StructureMatcher(ltol=ltol, stol=stol,
                                                                                         angle_tol=angle_tol).fit(
                    input_structure,
                    output_struct)
                return structure_matcher

            # filters through closest to base property in pymatgen if it's a value, based on previous list.
            if (type(base_property_py1) is float) or (type(base_property_py1) is int):
                for base_prop in base_property_reduced_and_hubbard:
                    if (type(base_prop) is int) or (type(base_prop) is float):
                        min_difference_rep.append(abs(base_prop - base_property_py1))
                        min_integers_rep.append(abs(base_prop - base_property_py1))
                    else:
                        min_difference_rep.append(base_prop)
                        continue
                # filters through structure matching and base property
                min_difference_rep_pos = [i for i, j in enumerate(min_difference_rep) if
                                          ((type(j) is float) or (type(j) is int))]
                min_difference_rep2 = []
                auid_reduced_and_hubbard2 = []
                for min_d in min_difference_rep_pos:
                    min_difference_rep2.append(min_difference_rep[min_d])
                    auid_reduced_and_hubbard2.append(auid_reduced_and_hubbard[min_d])

                assert len(min_difference_rep2) == len(auid_reduced_and_hubbard2)
                new_auid_reduced_and_hubbard = [x for _, x in
                                                sorted(zip(min_difference_rep2, auid_reduced_and_hubbard2),
                                                       key=lambda pair: pair[0])]
                # if no properties in AFLOW appear, disregards base property
                if len(new_auid_reduced_and_hubbard) == 0:
                    new_auid_reduced_and_hubbard = auid_reduced_and_hubbard
                if len(new_auid_reduced_and_hubbard) != 0:
                    for auid_reduced_and_h in new_auid_reduced_and_hubbard:
                        result_aflow_CONTCAR = aflow.search().filter(aflow.K.auid == auid_reduced_and_h)
                        for pA in result_aflow_CONTCAR:
                            species_CONT1 = pA.species
                            composition_CONT = pA.composition
                        substance_vector = []
                        species_CONT1[-1] = species_CONT1[-1].replace('\n', '')
                        for e in species_CONT1:
                            position_atom = species_CONT1.index(e)
                            substance_vector.append(e + str(composition_CONT[position_atom]))
                        perm_subs = permutations(substance_vector)
                        for spe in perm_subs:
                            spe_list = list(spe)
                            species_and_comp = ' '.join(spe_list)
                            species = ''.join([i for i in species_and_comp if not i.isdigit()])
                            res = []
                            for test_string in spe_list:
                                res.append(([str(i) for i in test_string if i.isdigit()]))
                            compo_list = [''.join(x) for x in res]
                            composition = ' '.join(compo_list)
                            is_structure_match = structurefiltering(species_CONT=species, comp_cont2=composition)
                            if is_structure_match is not np.bool_(True, dtype=bool):
                                continue
                            if is_structure_match is np.bool_(True, dtype=bool):
                                break
                        if os.path.exists('MATCOR_AFLOW_CONTCAR'):
                            os.remove("MATCOR_AFLOW_CONTCAR")
                        else:
                            print("MATCOR_AFLOW_CONTCAR not created")
                        if is_structure_match is not np.bool_(True, dtype=bool):
                            auid_AFLOW_base = 'NO_MATCH'
                            structrure_verification = 'NSV'
                            error_py_AFLOW = 'NSV'
                            continue
                        if is_structure_match is np.bool_(True, dtype=bool):
                            auid_AFLOW_base = auid_reduced_and_h
                            structrure_verification = 'SV'
                            break
                        break
                else:
                    auid_AFLOW_base = 'NO_MATCH'
                    error_py_AFLOW = 'NSV'
            # Value not provided by MP or base_property not a value:
            else:
                base_position_property = [i for i, j in enumerate(base_property_reduced_and_hubbard) if
                                          ((type(j) is float) or (type(j) is int))]
                base_property_and_hubbard2 = []
                auid_reduced_and_hubbard2 = []
                for min_d in base_position_property:
                    base_property_and_hubbard2.append(base_property_reduced_and_hubbard[min_d])
                    auid_reduced_and_hubbard2.append(auid_reduced_and_hubbard[min_d])

                assert len(base_property_and_hubbard2) == len(auid_reduced_and_hubbard2)
                new_auid_reduced_and_hubbard = [x for _, x in
                                                sorted(zip(base_property_and_hubbard2, auid_reduced_and_hubbard2),
                                                       key=lambda pair: pair[0])]
                if len(new_auid_reduced_and_hubbard) == 0:
                    new_auid_reduced_and_hubbard = auid_reduced_and_hubbard
                    base_property_and_hubbard2 = [0]
                if len(new_auid_reduced_and_hubbard) != 0:
                    for auid_reduced_and_h in new_auid_reduced_and_hubbard:
                        print(auid_reduced_and_h)
                        result_aflow_CONTCAR = aflow.search().filter(aflow.K.auid == auid_reduced_and_h)
                        for pA in result_aflow_CONTCAR:
                            species_CONT1 = pA.species
                            composition_CONT = pA.composition
                        substance_vector = []
                        species_CONT1[-1] = species_CONT1[-1].replace('\n', '')
                        for e in species_CONT1:
                            position_atom = species_CONT1.index(e)
                            substance_vector.append(e + str(composition_CONT[position_atom]))
                        perm_subs = permutations(substance_vector)
                        for spe in perm_subs:
                            spe_list = list(spe)
                            species_and_comp = ' '.join(spe_list)
                            species = ''.join([i for i in species_and_comp if not i.isdigit()])
                            res = []
                            for test_string in spe_list:
                                res.append(([str(i) for i in test_string if i.isdigit()]))
                            compo_list = [''.join(x) for x in res]
                            composition = ' '.join(compo_list)
                            is_structure_match = structurefiltering(species_CONT=species, comp_cont2=composition)
                            if is_structure_match is not np.bool_(True, dtype=bool):
                                continue
                            if is_structure_match is np.bool_(True, dtype=bool):
                                break
                        if os.path.exists('MATCOR_AFLOW_CONTCAR'):
                            os.remove("MATCOR_AFLOW_CONTCAR")
                        else:
                            print("MATCOR_AFLOW_CONTCAR not created")
                        if is_structure_match is not np.bool_(True, dtype=bool):
                            auid_AFLOW_base = 'NO_MATCH'
                            structrure_verification = 'NSV'
                            error_py_AFLOW = 'NSV'
                            continue
                        if is_structure_match is np.bool_(True, dtype=bool):
                            auid_AFLOW_base = auid_reduced_and_h
                            structrure_verification = 'SV'
                            break
                        break
                else:
                    auid_AFLOW_base = 'NO_MATCH'
                    error_py_AFLOW = 'NSV'
            # no match through any filtering
            if auid_reduced_and_sg == 'place_holder':
                auid_AFLOW_base = 'NO_MATCH'
                error_py_AFLOW = 'NO SG TO COMPARE'
            if len(auid_reduced_and_hubbard) == 0:
                auid_AFLOW_base = 'NO_MATCH'
                error_py_AFLOW = 'NO HUBBARD TO COMPARE'
            if len(auid_reduced_and_subs) == 0:
                auid_AFLOW_base = 'NO_MATCH'
                error_py_AFLOW = 'NO RESPONSE FROM AFLOW'
    if len(result_pymatgen) == 0:
        auid_AFLOW_base = 'NO_MATCH'
        error_py_AFLOW = 'NO RESPONSE FROM MP'
    return [auid_AFLOW_base, auid_reduced_and_hubbard, final_sg_aflow, space_group_pymatgen,
            structrure_verification, all_hubbards_aflow_reduced_and_hubbard, formula_input, error_py_AFLOW]


def AFLOW_to_pymatgen(auid_AFLOW, mpr, base_property_py1, base_property_py2, base_property_af, ltol, stol, angle_tol):
    # avoids possible reference errors
    key_to_pymatgen = ''
    final_sg_auid = ''
    material_id = 'place_holder'
    all_space_groups_py_final = []
    structure_verification = ''
    hubbard_AFLOW = ''
    mp_pool_hubbard_edited = []
    all_auid_sg0 = ''
    all_auid_sg = ''
    atoms = ''
    elements = ''
    prototype = ''
    formula_input = ''
    is_structure_match = ''
    base_property_AFLOW = ''
    species_CONT1 = ''
    material_ids_hubbard = []
    all_space_groups_hubbard = []
    all_base_prop_hubbard = []
    all_space_groupsPY = []
    material_ids = []
    all_icsd_ids = []
    substances_AFLOW = []
    min_integers = []
    all_base_prop = []
    all_space_group_pymatgen = []
    is_hubbard_py = []
    match_representative = []
    match_base = []
    min_integers = []
    all_space_groups_representative = []
    mp_pool_hubbard = []
    composition_CONT = ''
    min_difference_representative = []
    match_structure = []
    error_aflow_py = ''
    try:
        result_aflow = aflow.search().filter(aflow.K.auid == auid_AFLOW)
        for entry in result_aflow:
            prototype = entry.prototype
            prototype = prototype.replace('0.', '')
            atoms = entry.composition
            elements = entry.species
            all_auid_sg = entry.sg2
            all_auid_sg0 = entry.sg
            base_property_AFLOW = getattr(entry, base_property_af)
            hubbard_AFLOW = entry.ldau_TLUJ
            formula_input = str(entry.compound)
        # sg for after calculation
        auid_sg_after0 = all_auid_sg0[2]
        auid_sg_after0 = auid_sg_after0.replace('{', '')
        auid_sg_after0 = auid_sg_after0.replace('}', '')
        auid_sg_after0 = auid_sg_after0.replace('\n', '')

        auid_sg_after = all_auid_sg[2]
        auid_sg_after = auid_sg_after.replace('{', '')
        auid_sg_after = auid_sg_after.replace('}', '')
        auid_sg_after = auid_sg_after.replace('\n', '')

        # sg2 for before calculation
        auid_sg_before = all_auid_sg[0]
        auid_sg_before = auid_sg_before.replace('{', '')
        auid_sg_before = auid_sg_before.replace('}', '')
        auid_sg_before = auid_sg_before.replace('\n', '')

        auid_sg_before0 = all_auid_sg0[0]
        auid_sg_before0 = auid_sg_before0.replace('{', '')
        auid_sg_before0 = auid_sg_before0.replace('}', '')
        auid_sg_before0 = auid_sg_before0.replace('\n', '')

        final_sg_auid = str(auid_sg_before) + ',' + str(auid_sg_after) + ',' + str(auid_sg_before0) + ',' + str(
            auid_sg_after0)
        auid_sg_after = [auid_sg_before, auid_sg_after, auid_sg_before0, auid_sg_after0]
        # changes AFLOW to pymatgen format
        result = atoms[0]
        for i in atoms[1:]:
            result = gcd(result, i)
        pretty_natoms = []
        for x in atoms:
            pretty_natoms.append(int(x / result))
        for n, s in enumerate(pretty_natoms):
            if s == 1:
                pretty_natoms[n] = ' dummy'
        elements[-1] = elements[-1].replace('\n', '')
        substance_vector = []
        substance_vector1 = []
        for e in elements:
            position_atom = elements.index(e)
            substance_vector.append(e + str(atoms[position_atom]))
            substance_vector1.append(e + str(pretty_natoms[position_atom]))
        substance_full = ''.join(substance_vector1)
        if 'dummy' in substance_full:
            substance_full = substance_full.replace(' dummy', '1')
        prototype = prototype[prototype.find('ICSD'):]
        prototype = prototype.replace('ICSD_', '')
        icsd_AFLOW_vector = str(prototype.replace(substance_full, ''))
        icsd_AFLOW = str(''.join(i for i in icsd_AFLOW_vector if i.isdigit()))
        icsd_AFLOW = icsd_AFLOW.replace('\n', '')

        perm = permutations(substance_vector1)
        for i in perm:
            j = str(i)
            j = j.replace(' dummy', '')
            j = j.replace('(', '')
            j = j.replace("'", '')
            j = j.replace(')', '')
            j = j.replace(', ', '')
            substances_AFLOW.append(j)
        # finds all possible icsds ids
    except:
        material_id = 'NO_MATCH'
        error_aflow_py = 'NO RESPONSE FROM AFLOW'
    for substance in substances_AFLOW:
        substance = substance.replace(',', '')
        criteria1 = {'pretty_formula': substance}
        properties1 = supported_properties_py
        result_pymatgen = mpr.query(criteria=criteria1, properties=properties1)
        if len(result_pymatgen) != 0:
            key_to_pymatgen = 'Open'
            formula_input = substance
            for elem in result_pymatgen:
                property_base = elem.get(base_property_py1)
                all_icsd_ids.append(elem.get('icsd_ids'))
                material_ids.append(elem.get('material_id'))
                all_space_groupsPY.append(elem.get('spacegroup'))
                is_hubbard_py.append(elem.get('is_hubbard'))

                if len(base_property_py2) == 0:
                    all_base_prop.append(elem.get(base_property_py1))
                if len(base_property_py2) != 0:
                    if property_base is None:
                        all_base_prop.append('None')
                    else:
                        all_base_prop.append(elem.get(base_property_py1).get(base_property_py2))
            for space_group in all_space_groupsPY:
                space_group_pymatgen = str(space_group.get('symbol')) + ' #' + str(space_group.get('number'))
                all_space_group_pymatgen.append(space_group_pymatgen)
            break
        if len(result_pymatgen) == 0:
            key_to_pymatgen = 'Closed'
            continue
    if key_to_pymatgen == 'Open':
        # verifies MP returns ICSD
        lenght = []
        for i in all_icsd_ids:
            lenght.append(len(i))
        max_icsd = max(lenght)
        if hubbard_AFLOW is None:
            all_reduced_and_hubbard_pos_true = [i for i, j in enumerate(is_hubbard_py) if
                                                j is False]
            if len(all_reduced_and_hubbard_pos_true) != 0:
                for hubbard_pos in all_reduced_and_hubbard_pos_true:
                    material_ids_hubbard.append(material_ids[hubbard_pos])
                    all_space_groups_hubbard.append(all_space_group_pymatgen[hubbard_pos])
                    all_base_prop_hubbard.append(all_base_prop[hubbard_pos])
            if len(all_reduced_and_hubbard_pos_true) == 0:
                material_ids_hubbard = material_ids
                all_space_groups_hubbard = all_space_group_pymatgen
                all_base_prop_hubbard = all_base_prop
        if hubbard_AFLOW is not None:
            all_reduced_and_hubbard_pos_false = [i for i, j in enumerate(is_hubbard_py) if
                                                 j is True]
            if len(all_reduced_and_hubbard_pos_false) != 0:
                for hubbard_neg in all_reduced_and_hubbard_pos_false:
                    material_ids_hubbard.append(material_ids[hubbard_neg])
                    all_space_groups_hubbard.append(all_space_group_pymatgen[hubbard_neg])
                    all_base_prop_hubbard.append(all_base_prop[hubbard_neg])
            if len(all_reduced_and_hubbard_pos_false) == 0:
                material_ids_hubbard = material_ids
                all_space_groups_hubbard = all_space_group_pymatgen
                all_base_prop_hubbard = all_base_prop
        # for G_VRH error: prioritizes matching of integers and not Hubbard U
        check_base_result = [s for s in all_base_prop_hubbard if s is not None]
        if len(check_base_result) == 0:
            material_ids_hubbard = material_ids
            all_space_groups_hubbard = all_space_group_pymatgen
            all_base_prop_hubbard = all_base_prop
        if len(all_space_groups_hubbard) != 0:
            def representative():
                all_position = []
                for sg_py in all_space_groups_hubbard:
                    if sg_py in auid_sg_after:
                        all_position = [i for i, j in enumerate(all_space_groups_hubbard) if
                                        j in sg_py]
                if len(all_position) != 0:
                    for matchingpos in all_position:
                        match_representative.append(material_ids_hubbard[matchingpos])
                        match_base.append(all_base_prop_hubbard[matchingpos])
                        all_space_groups_representative.append(all_space_groups_hubbard[matchingpos])
                        # further filters through base property, if value
                    if (type(base_property_AFLOW) is float) or (type(base_property_AFLOW) is int):
                        for base_proper in match_base:
                            if (type(base_proper) is int) or (type(base_proper) is float):
                                min_difference_representative.append(abs(base_proper - base_property_AFLOW))
                                min_integers.append(abs(base_proper - base_property_AFLOW))
                                min_difference_val = min(min_integers)
                            else:
                                min_difference_representative.append(base_proper)
                                continue

                all_space_groups_py_final = all_space_groups_representative
                print('Entries with same sg and formula, and potential hubbard U selection')
                print(match_representative)
                return [match_representative, all_space_groups_py_final, min_difference_representative]

            result_representative = representative()
            match_structure = result_representative[0]
            all_space_groups_py_final = result_representative[1]
            min_difference = result_representative[2]

            # gets hubbard for log file
            for pool_mp in match_structure:
                criteria_hubbard = {'material_id': pool_mp}
                properties_hubbard = supported_properties_py  # do not change
                result_pymatgen_hubbard = mpr.query(criteria=criteria_hubbard, properties=properties_hubbard)
                for hubres in result_pymatgen_hubbard:
                    hubbard_val_py = str(hubres.get('hubbards'))
                    hubbard_val_py = hubbard_val_py.replace('\n', '')
                    mp_pool_hubbard.append(hubbard_val_py)
            mp_pool_hubbard_edited = []
            for mpedit in mp_pool_hubbard:
                if (len(mpedit) == 0) or (mpedit == '{}'):
                    mp_pool_hubbard_edited.append('None')
                else:
                    mp_pool_hubbard_edited.append(mpedit)
            # filters through structure matching; if there is no base property available, it will ignore this
            min_difference_rep_pos = [i for i, j in enumerate(min_difference) if
                                      ((type(j) is float) or (type(j) is int))]
            min_difference_rep2 = []
            match_structure2 = []
            for min_d in min_difference_rep_pos:
                min_difference_rep2.append(min_difference[min_d])
                match_structure2.append(match_structure[min_d])
            assert len(min_difference_rep2) == len(match_structure2)
            new_match_structure = [x for _, x in
                                   sorted(zip(min_difference_rep2, match_structure2), key=lambda pair: pair[0])]
            if len(new_match_structure) == 0:
                new_match_structure = match_structure

            def CONTCAR_to_pymatgen(species_CONT, comp_cont2):
                a0 = result_aflow_CONTCAR[0]
                files = a0.files
                if 'CONTCAR.relax' in files:
                    pos_file = files.index('CONTCAR.relax')
                    MATCOR_tar = files[pos_file]('MATCOR_AFLOW_CONTCAR')
                else:
                    print('CONTCAR.relax file not provided by AFLOW')
                # makes changes for CONTCAR file in MP format
                with open('MATCOR_AFLOW_CONTCAR', 'r') as f:
                    text = f.readlines()
                    text[5] = species_CONT + '\n' + comp_cont2 + '\n'
                    if 'Direct' not in str(text[6]):
                        text.pop(7)
                my_file = open('MATCOR_AFLOW_CONTCAR', "w")
                new_file_contents = "".join(text)
                my_file.write(new_file_contents)
                my_file.close()
                output_struct = Structure.from_file('MATCOR_AFLOW_CONTCAR')
                return output_struct

            # gets CONTCAR file in AFLOW ready for structure matching in pymatgen format
            result_aflow_CONTCAR = aflow.search().filter(aflow.K.auid == auid_AFLOW)
            if len(new_match_structure) != 0:
                for mp_id in new_match_structure:
                    output_struct = mpr.get_structure_by_material_id(mp_id)
                    for pA in result_aflow_CONTCAR:
                        species_CONT1 = pA.species
                        composition_CONT = pA.composition
                    substance_vector = []
                    species_CONT1[-1] = species_CONT1[-1].replace('\n', '')
                    for e in species_CONT1:
                        position_atom = species_CONT1.index(e)
                        substance_vector.append(e + str(composition_CONT[position_atom]))
                    perm_subs = permutations(substance_vector)
                    for spe in perm_subs:
                        spe_list = list(spe)
                        species_and_comp = ' '.join(spe_list)
                        species = ''.join([i for i in species_and_comp if not i.isdigit()])
                        res = []
                        for test_string in spe_list:
                            res.append(([str(i) for i in test_string if i.isdigit()]))
                        compo_list = [''.join(x) for x in res]
                        composition = ' '.join(compo_list)
                        input_structure = CONTCAR_to_pymatgen(species_CONT=species, comp_cont2=composition)
                        is_structure_match = pymatgen.analysis.structure_matcher.StructureMatcher(ltol=ltol, stol=stol,
                                                                                                  angle_tol=angle_tol).fit(
                            input_structure,
                            output_struct)
                        if is_structure_match is not np.bool_(True, dtype=bool):
                            continue
                        if is_structure_match is np.bool_(True, dtype=bool):
                            break
                    if is_structure_match is not np.bool_(True, dtype=bool):
                        material_id = 'NO_MATCH'
                        structure_verification = 'NSV'
                        error_aflow_py = 'NSV'
                        continue
                    if is_structure_match is np.bool_(True, dtype=bool):
                        material_id = mp_id
                        structure_verification = 'SV'
                        break
            else:
                material_id = 'NO_MATCH'
                error_aflow_py = 'NSV'
        if len(all_space_groups_hubbard) == 0:
            material_id = 'NO_MATCH'
            error_aflow_py = 'NO SG TO COMPARE'
    else:
        material_id = 'NO_MATCH'
        error_aflow_py = 'NO RESPONSE FROM MP'
    assert len(match_structure) == len(all_space_groups_py_final) == len(mp_pool_hubbard_edited)

    return [material_id, match_structure, all_space_groups_py_final, final_sg_auid, structure_verification,
            mp_pool_hubbard_edited, formula_input, error_aflow_py]


if os.path.isfile('./' + MATCOR_out) is True:
    file_okay = str(input('The ' + MATCOR_out + ' file already exists. Do you wish to continue (Y/N)? '))
    print(os.path.dirname(os.path.abspath(MATCOR_out)))
    if file_okay == 'Y':
        key_MATCOR = True
    if file_okay == 'N':
        key_MATCOR = False
else:
    key_MATCOR = True
minimum_difference_list = []

minimum_difference = ''
test_mp = ''
if key_MATCOR is True:
    with open(MATCOR_out, 'w') as f:
        f.write('MATCOR: A Program for the Cross-Validation of Material Properties Between Databases\n')
        f.write('\nSV == Structure VERIFIED with pymatgen structure matcher\n')
        f.write('NSV == Structure NOT VERIFIED with pymatgen structure matcher\n')
        f.write('HUB == Hubbard-U used\n')
        f.write('DFT == Hubbard-U not used\n' + '\n')
        run = 0
        assert (type(base_property_MP_1) is not int) or (type(base_property_MP_2) is not int) or (
                type(base_property_AFLOW) is not int), 'Select a property supported by ' \
                                                       'Materials Project and AFLOW '
        assert (len(base_property_MP_1) != 0) and (len(base_property_AFLOW) != 0), 'Select a property supported by ' \
                                                                                   'Materials Project and AFLOW '
        assert base_property_MP_1 in supported_properties_py, 'Property selected must be a supported property by ' \
                                                              'Materials Project '

        assert len(materials) != 0, 'List of materials must contain at least 1 material'
        blank_space_output = ''
        for id in materials:
            id.replace('\n', '')
            assert ('mp-' in id) or ('mvc-' in id) or ('aflow:' in id), 'Format of id is incorrect'
            key_headers = key_headers + 1
            try:
                if ('mp-' in id) or ('mvc' in id):
                    result1 = pymatgen_to_AFLOW(id=id, mpr=mpr, base_property_py1=base_property_MP_1,
                                                base_property_py2=base_property_MP_2,
                                                base_property_af=base_property_AFLOW, ltol=ltol_in, stol=stol_in,
                                                angle_tol=angle_tol_in)
                    material_id = id
                    auid_AFLOW_base = result1[0]
                    formula_input = result1[6]
                    pool_possible_matches = result1[1]
                    struct_verification = result1[4]
                    output_material = auid_AFLOW_base
                    output_error = result1[7]

                if 'aflow' in id:
                    result2 = AFLOW_to_pymatgen(auid_AFLOW=id, mpr=mpr, base_property_py1=base_property_MP_1,
                                                base_property_py2=base_property_MP_2,
                                                base_property_af=base_property_AFLOW, ltol=ltol_in, stol=stol_in,
                                                angle_tol=angle_tol_in)
                    auid_AFLOW_base = id
                    material_id = result2[0]
                    formula_input = result2[6]
                    struct_verification = result2[4]
                    pool_possible_matches = result2[1]
                    output_material = material_id + ((22 - len(
                        material_id)) * ' ')
                    output_error = result2[7]
                if (auid_AFLOW_base == 'NO_MATCH') or (material_id == 'NO_MATCH') or (len(auid_AFLOW_base) == 0) or (
                        len(material_id) == 0) or (len(pool_possible_matches) == 0):
                    best_error = 'pool_best_' + str(materials.index(id)) + ':'
                    f.write(best_error + (30 - len(best_error)) * ' ')
                    f.write(id)
                    best_error_id = (36 - len(best_error + id)) * ' '
                    f.write(best_error_id)
                    if output_error == '':
                        f.write('NO_MATCH \n')
                    else:
                        f.write('NO_MATCH   ' + output_error + '\n')
                    continue
                if (auid_AFLOW_base != 'NO_MATCH') and (material_id != 'NO_MATCH') and (
                        len(pool_possible_matches) != 0):
                    result_aflow2 = aflow.search().filter(aflow.K.auid == auid_AFLOW_base)
                    for entry2 in result_aflow2:
                        property_AFLOW_email_0 = getattr(entry2, base_property_AFLOW)
                        if (type(property_AFLOW_email_0) is int) or (type(property_AFLOW_email_0) is float):
                            property_AFLOW_email = round(Decimal(property_AFLOW_email_0), 4)
                        else:
                            property_AFLOW_email = property_AFLOW_email_0
                        property_AFLOW = str(property_AFLOW_email)
                        hubbard_val_aflow = str(entry2.ldau_TLUJ)
                    if hubbard_val_aflow == 'None':
                        hubbard_U_aflow = 'DFT'
                        hubbard_f_aflow = 'None'
                    else:
                        hubbard_U_aflow = 'HUB'
                        hubbard_f_aflow = hubbard_val_aflow
                    criteria1 = {'material_id': material_id}
                    properties1 = [base_property_MP_1, 'pretty_formula', 'is_hubbard', 'hubbards']
                    materials_py = mpr.query(criteria=criteria1, properties=properties1)
                    hubbard_val_py = ''
                    for material in materials_py:
                        property_mp = material.get(base_property_MP_1)
                        hubbard_f_py = str(material.get('is_hubbard'))
                        if hubbard_f_py == 'False':
                            hubbard_U_py = 'DFT'
                        else:
                            hubbard_U_py = 'HUB'
                        hubbard_f_py = str(hubbard_f_py.replace('\n', ''))
                        hubbard_val_py = str(material.get('hubbards'))
                        hubbard_val_py = hubbard_val_py.replace('\n', '')
                        if len(base_property_MP_2) == 0:
                            if (type(property_mp) is int) or (type(property_mp) is float):
                                property_py_email = round(Decimal(property_mp), 4)
                                property_py = str(property_py_email)
                            else:
                                print('yoy got it')
                                property_py_email = property_mp
                                property_py = property_py_email
                            if property_mp is None:
                                property_py = 'None'
                                property_py_email = 'None'

                        if len(base_property_MP_2) != 0:
                            test_mp = material.get(base_property_MP_1)
                            if test_mp is not None:
                                property_mp = material.get(base_property_MP_1).get(base_property_MP_2)
                            if test_mp is None:
                                property_mp = material.get(base_property_MP_1)
                            if (type(property_mp) is int) or (type(property_mp) is float):
                                property_py_email = round(Decimal(property_mp), 4)
                                property_py = str(property_py_email)
                            else:
                                property_py_email = property_mp
                                property_py = property_py_email
                            if property_mp is None:
                                property_py = 'None'
                                property_py_email = 'None'
                    if (type(property_py_email) is float or type(property_py_email) is int or type(
                            property_py_email) is decimal.Decimal) and (
                            type(property_AFLOW_email) is float or type(property_AFLOW_email) is int or type(
                        property_AFLOW_email) is decimal.Decimal):
                        minimum_difference = abs(round(Decimal(float(property_mp) - float(property_AFLOW_email)), 4))
                    else:
                        minimum_difference = 'not a value'
                    minimum_difference_list.append(minimum_difference)
                    mydata = [
                        [material_id, str(property_py), str(hubbard_f_py), auid_AFLOW_base, property_AFLOW,
                         hubbard_f_aflow]]
                    y = ''
                    headers = ''
                    if key_headers == 1:
                        headers = ['material_id', base_property_MP_1, 'hubbard_py', 'auid_AFLOW_base',
                                   base_property_AFLOW, 'hubbard_aflow']
                    if key_headers != 1:
                        headers = [y.ljust(len('material_id')), y.ljust(len(base_property_MP_1)),
                                   y.ljust(len('hubbard_py')),
                                   y.ljust(len('auid_AFLOW_base')), y.ljust(len(base_property_AFLOW)),
                                   y.ljust(len('hubbard_aflow'))]
                    print(tabulate(mydata, headers=headers))
                    # creates log file
                    run = run + 1
                    pool_hubbards = []
                    if len(pool_possible_matches) != 0:
                        if ('mp-' in id) or ('mvc' in id):
                            hubbard_U = hubbard_U_py + '-' + hubbard_U_aflow
                            if hubbard_val_py == '{}':
                                input_hubbard = 'None'
                            else:
                                input_hubbard = hubbard_val_py
                            all_sg_pool = result1[2]
                            sg_input = result1[3]
                            struct_verification = result1[4]
                            pool_hubbards = result1[5]
                            blank_space_input = (25 - len(id)) * ' '
                            blank_space_pool = 3 * ' '

                        if 'aflow' in id:
                            hubbard_U = hubbard_U_aflow + '-' + hubbard_U_py
                            input_hubbard = hubbard_val_aflow
                            all_sg_pool = result2[2]
                            sg_input = result2[3]
                            blank_space_input = 3 * ' '
                            struct_verification = result2[4]
                            pool_hubbards = result2[5]  # missing for aflow

                            blank_space_pool = 'length of matching id'
                        header_id = 'pool_' + str(materials.index(id)) + '_' + '0: '
                        f.write(str(
                            header_id + id) + blank_space_input + sg_input + ((70 - len(
                            sg_input)) * ' ') + input_hubbard + blank_space_input + formula_input + '\n')
                        blank_space_output = str(header_id + id + blank_space_input)
                        run_in = 1
                        sg_run = 0
                        if len(pool_possible_matches) != 0:
                            for match in pool_possible_matches:
                                if run_in <= len(pool_possible_matches):
                                    run_in = run_in + 1
                                    matches_pool = 'pool_' + str(materials.index(id)) + '_' + str(
                                        (pool_possible_matches.index(match)) + 1)
                                    if blank_space_pool == 'length of matching id':
                                        f.write(matches_pool + ': ' + match + (
                                                len(blank_space_output) - len(matches_pool + ': ' + match)) * ' ')
                                        f.write(
                                            str(all_sg_pool[sg_run]) + ((len(str(
                                                header_id + id) + blank_space_input + sg_input + ((70 - len(
                                                sg_input)) * ' ')) - len((matches_pool + ': ' + match + (
                                                    len(blank_space_output) - len(
                                                matches_pool + ': ' + match)) * ' ') + str(
                                                all_sg_pool[sg_run]))) * ' ') + str(
                                                pool_hubbards[sg_run]) + '\n')
                                        sg_run = sg_run + 1
                                    else:
                                        f.write(matches_pool + ': ' + match + (
                                                len(blank_space_output) - len(matches_pool + ': ' + match)) * ' ')
                                        f.write(
                                            str(all_sg_pool[sg_run]) + ((len(str(
                                                header_id + id) + blank_space_input + sg_input + ((70 - len(
                                                sg_input)) * ' ')) - len((matches_pool + ': ' + match + (
                                                    len(blank_space_output) - len(
                                                matches_pool + ': ' + match)) * ' ') + str(
                                                all_sg_pool[sg_run]))) * ' ') + str(
                                                pool_hubbards[sg_run]) + '\n')
                                        sg_run = sg_run + 1
                            hubbard_val_aflow = '"' + hubbard_val_aflow + '"'
                            labelLine = list()
                            valueLine = list()
                            if ('mp-' in id) or ('mvc' in id):
                                material_id = 'pool_best_' + str(materials.index(id)) + ':' + ' ' + str(
                                    formula_input) + str(
                                    ((22 - len('pool_best_' + str(materials.index(id)) + ':' + ' ' + str(
                                        formula_input))) * ' ')) + str(id) + str(((22 - len(id)) * ' '))
                            if 'aflow' in id:
                                material_id = 'pool_best_' + str(materials.index(id)) + ':' + ' ' + str(
                                    formula_input) + str(
                                    ((22 - len('pool_best_' + str(materials.index(id)) + ':' + ' ' + str(
                                        formula_input))) * ' ')) + str(id) + str(((22 - len(id)) * ' '))

                            data = [['input_id' + str(((30 - len(material_id)) * ' ')),
                                     'hubbard_U' + str(((15 - len(material_id)) * ' ')),
                                     'SV//NSV' + str(((15 - len(material_id)) * ' ')),
                                     'property_py' + str(((15 - len(material_id)) * ' ')),
                                     'output_material' + str(((15 - len(material_id)) * ' ')),
                                     'property_AFLOW' + str(((15 - len(material_id)) * ' ')),
                                     'min_difference'],
                                    [material_id, str(hubbard_U), str(struct_verification), str(property_py),
                                     str(output_material),
                                     str(property_AFLOW), (str(minimum_difference)).replace('\n', '')]]
                            for label, value in zip(
                                    *data):  # unzips the list, so each elements (label and value) get extracted pairwise
                                padding = max(len(str(label)),
                                              len(str(value)))  # what is longer, the label or the value?
                                labelLine.append(
                                    '{0:<{1}}'.format(label,
                                                      padding))  # generate a string with the variable whitespace padding
                                valueLine.append(
                                    '{0:<{1}}'.format(value,
                                                      padding))  # generate a string with the variable whitespace padding
                            labelLine = [label for label in labelLine]
                            valueLine = [value for value in valueLine]
                            f.write('\t'.join(labelLine) + '\n')
                            f.write('\t'.join(valueLine) + '\n')
            except:
                print(id)
                print('MATCOR did not respond for this entry.')
                f.write('NO_MATCH        MATCOR ERROR \n')

    # sends email
    EMAIL_ADRESS = 'your adress'
    EMAIL_PASSWORD = 'your password'
    new_minimum_difference_pos = [i for i, j in enumerate(minimum_difference_list) if
                                      abs(j) > max_difference]
    # uncomment here to send email
    #     if len(new_minimum_difference_pos) != 0:
    #         mail_content = 'MATCOR is an ' \
    #                                'open-source user-friendly, easily adaptable software to facilitate the data curation process in entries ' \
    #                                'between AFLOW and Materials Project. By the choice of the user, an automated email was generated to ' \
    #                                'notify both MP and AFLOW about specific entries whose difference value for a specific property were ' \
    #                                'greater than a user-defined difference.\n\n\n\n'
    #         mail_content = mail_content + 'Attached is the output file produced by MATCOR'
    #         # The mail addresses and password
    #         sender_address = 'your address'
    #         sender_pass = 'your password'
    #         receiver_address = 'feedback@materialsproject.org', or 'aflow@groups.io'
    #         # Setup the MIME
    #         message = MIMEMultipart()
    #         message['From'] = sender_address
    #         message['To'] = receiver_address
    #         message['Subject'] = 'MATCOR: Significant differences found between AFLOW and Materials Project'
    #         # The subject line
    #         # The body and the attachments for the mail
    #         message.attach(MIMEText(mail_content, 'plain'))
    #         attach_file_name = MATCOR_out
    #         attach_file = open(attach_file_name, 'rb')  # Open the file as binary mode
    #         payload = MIMEBase('application', 'octate-stream')
    #         payload.set_payload((attach_file).read())
    #         encoders.encode_base64(payload)  # encode the attachment
    #         # add payload header with filename
    #         payload.add_header('MATCOR', 'attachment', filename=attach_file_name)
    #         message.attach(payload)
    #         # Create SMTP session for sending the mail
    #         session = smtplib.SMTP('smtp.gmail.com', 587)  # use gmail with port
    #         session.starttls()  # enable security
    #         session.login(sender_address, sender_pass)  # login with mail_id and password
    #         text = message.as_string()
    #         session.sendmail(sender_address, receiver_address, text)
    #         session.quit()
    #         print('Mail Sent')
    #     else:
    #             print('Properties not values, not found, or within established difference. Email not sent')
else:
    print('Create a new txt file')
