import re


class simReader:
    """
    this class reads file to family list of IDs and list of alleles
    file columns: family ID, seq in family (F/M/child num), 4 cols for each gene: allele-1*, allele-2*,allele-1,allele-2
    """
    def __init__(self, data_path):
        self._f = open(data_path, 'rt')
        dummy = next(self._f)
        self._curr_line = (next(self._f)).split(',')
        self._currIndex = 1
        self._end = False

    @staticmethod
    def extract_als(al_list, amb, alllele):
        """
        extract als in format to work with, ignoring letters that represent ambiguity.
        first 2 are updated and will be taken if exist.
        :param al_list: from file
        :return: the 2 alleles that were extracted as the best of 4 in file.
        """
        al1 = al_list[0]
        al2 = al_list[1]
        if al_list[0] == '':  # no high res
            al1 = al_list[2]
            al2 = al_list[3]
        if al2 == '':  # homozygous
            al2 = al1
        al1 = str(al1).replace('p', '')
        al2 = str(al2).replace('p', '')
        if ":" in al1 and al1.split(':')[1].isupper():
            amb[alllele + '*' + al1.split(':')[0]] = al1.split(':')[1]  # {'A*02': 'BJFVY'}
        if ":" in al2 and al2.split(':')[1].isupper():
            amb[alllele + '*' + al2.split(':')[0]] = al2.split(':')[1]
        new_al1 = re.sub('[a-zA-Z]', '', al1)  # remove ambiguity
        new_al2 = re.sub('[a-zA-Z]', '', al2)  # remove ambiguity
        if new_al1 == '\n':
            new_al1 = ''
        if new_al2 == '\n':
            new_al2 = ''
        return new_al1, new_al2

    @staticmethod
    def parents_num(family_ids):
        """
        :param family_ids: ids list of family
        :return: num of parents
        """
        if 'M' in family_ids and 'F' in family_ids:
            return 2
        elif 'M' in family_ids or 'F' in family_ids:
            return 1
        return 0

    def get_family(self):
        """
        read from current place in file all rows with the family id
        :return: list of ids and all lists of alleles for each
        """
        if self._end:
            print('end of file')
            return None, None, None, None
        fam_id = self._curr_line[0]
        i = self._currIndex
        family_ids = [fam_id]
        family_als = []
        amb = {}

        seq, ind = self.get_individual(self._curr_line, amb)
        family_als.append(ind)
        family_ids.append(seq)
        self._end = True
        for line in self._f:
            self._end = False
            if line == '\n':
                self._end = True
                break
            line_list = line.split(',')
            if line_list[0] == fam_id:
                seq, ind = self.get_individual(line_list, amb)
                if seq != "empty_member":  # if its not an empty list --> change to None?
                    # identical child - ignore from one of them
                    index_ind = family_als.index(ind) if ind in family_als else -1
                    if index_ind == -1 or not (isinstance(seq, (int, float))) \
                            or not (isinstance(family_ids[index_ind + 1], (int, float))):
                        family_als.append(ind)
                        family_ids.append(seq)
            else:
                self._curr_line = line_list
                par_num = self.parents_num(family_ids)
                return family_ids, family_als, par_num, amb
        if self._end:
            return None, None, None, None
        par_num = self.parents_num(family_ids)
        return family_ids, family_als, par_num, amb

    def get_individual(self, line, amb):
        """
        create list of alleles for individual from row in sheet.
        :param line: of individual's data.
        :return: individuals id and list of alleles.
        """
        seq = line[1]  # individual seq number ID
        ind = []
        if len(line[-1]) > 1 and line[-1][-1] == '\n':
            line[-1] = line[-1][:-1]
        temp_d = {2: 'A', 6: 'B', 10: 'C', 14: 'DRB1', 18: 'DQB1'}
        for i in range(2, 21, 4):
            val = line[i:i+4]
            data_pair = list(self.extract_als(val, amb, temp_d[i]))
            ind.append(data_pair)  # adding individual's alleles
        if not any(item for sublist in ind for item in sublist):  # completely empty: [("",""),("","")..]
            seq = None  # no data about this member
        return seq, ind
