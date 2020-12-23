from GR_code.GG_GRAMM.code.als import Als


def is_valid(fam_d, al_types, par_num):
    """
    validation tests of alleles data
    return: 0 if valid. string with error description otherwise.
    """
    for types in al_types:
        for key in fam_d:
            for single_al in fam_d[key][types]:
                if single_al == "" or single_al == " ":
                    continue
                elif ":" in single_al:
                    parts = single_al.split(":")
                    for part in parts:
                        if not (part.isnumeric() or part.isupper() or part == ''):
                            return "there is an allele which contains invalid character " \
                                   "(must contains numbers or uppercase only)"
                else:
                    if not (single_al.isnumeric() or single_al.isupper()):
                        return "there is an allele which contains invalid character " \
                               "(must contains numbers or uppercase only)"
    if len(fam_d) < 2:
        return "missing data for family (less than 2 people)"
    if len(fam_d) == 2 and par_num == 2:
        return "missing data for family - no children"
    for types in al_types:
        lst = Als()
        for key in fam_d:  # F, M, 1, 2...
            if any(fam_d[key][types]):  # not empty
                lst = fam_d[key][types].merge(lst)
        if len(lst) > 4:
            return "too many alleles in family"
    if par_num == 2:
        for types in al_types:
            fm_als = fam_d['F'][types] + fam_d['M'][types]
            for key in fam_d:
                if key != 'F' and key != 'M' and len(fm_als) == 4 and all(fm_als):
                    in_fm = fam_d[key][types].sub_lst(fm_als)
                    if not in_fm:
                        return "there is an allele in child that doesn't exist in parents"

    return 0


def list_to_dict(family_ids, family_als, als_names):
    """
    convert family data from list to dict
    :param family_ids: types of family members (F, M, 1, 2...)
    :param family_als: alleles data of family
    :param als_names: alleles types (A, B, C, DR, DQ)
    :return: dict with data about the family
    """
    fam_dict = {}
    for i in range(len(family_ids)):
        fam_member = family_ids[i]
        fam_dict[fam_member] = {}
        for j in range(len(family_als[0])):
            als_data = Als()
            als_data.append(family_als[i][j][0])
            als_data.append(family_als[i][j][1])
            fam_dict[fam_member][als_names[j]] = als_data
    return fam_dict


def child_ls_merge(child_ls, child_d, types):
    """
    merge alleles of children to list (ex. [01:03, 08, 01, 08] -> [01:03, 08])
    @param child_ls: child list
    @param child_d: child dictionary
    @param types: A/B/C/DR/DB
    @return: merged list
    """
    lst = Als()
    for child in child_ls:
        if any(child_d[child][types]):
            lst = child_d[child][types].merge(lst)
    return lst


# convert parents chromosomes to list (for writing in file)
def options_num_fm(father, mother):
    op_lst = []
    i = 0
    for chro in [father.ch1, father.ch2, mother.ch1, mother.ch2]:
        cur_op = 1
        for value in chro.values():
            if len(value) > 1:
                cur_op = cur_op * len(value)
        op_lst.append(cur_op)
        i += 1
    return op_lst


# write the current error to file
def errors_report(f_er, count, ind, error_name, fam_id):
    count[ind] += 1
    f_er.write('family "' + str(fam_id) + '": ' + error_name + '\n')  # write error to errors file


def sum_errors(f, er_count, fam_count, inconcis_er, all_err):
    f.write('\n\n***** Errors Summary *****' + '\n')
    f.write('Validation errors in pre-processing: ' + str(er_count[0]) +
            ' families. ' + "{:.1f}".format(er_count[0] / fam_count * 100) + '%' + ' from data\n')
    f.write('Can not embed child in parents data: ' + str(er_count[1]) +
            ' families. ' + "{:.1f}".format(er_count[1] / fam_count * 100) + '%' + ' from data\n')
    f.write('Errors in adding data. contradiction between parents and child: ' + str(er_count[2]) +
            ' families. ' + "{:.1f}".format(er_count[2] / fam_count * 100) + '%' + ' from data\n')
    f.write('Errors in creating GL string: ' + str(er_count[3]) +
            ' parents. ' + "{:.1f}".format(er_count[3] / fam_count * 50) + '%' + ' from data\n')
    f.write('Errors in inconsistency between parents and children: ' + str(inconcis_er) + ' families. ' +
            "{:.1f}".format(inconcis_er / int(fam_count) * 100) + '% from data\n')
    f.write('\nTotal errors:\n' + str(all_err) + ' from ' + str(fam_count) + ' families\n')
