__author__ = 'walzer'
__copyright__ = "M. Walzer"
__license__ = "BSD"
__maintainer__ = "walzer"
__email__ = "walzer<at>informatik.uni-tuebingen.de"

import xlrd
import xlwt
from xlutils.copy import copy
import string


def xlsDictReader(f, sheet_index=0, sheet_name='', no_header=False):
    """
    reads a sheet of a xls file into a list of dicts
    :param f: string for the filename
    :param sheet_index: the index of the sheet to access
    :param sheet_name: alternative to the sheet_index if the sheets in f are named
    :param no_header: set to true if the targeted sheet has no headers
    """
    book = xlrd.open_workbook(f)
    sheet = None
    if sheet_name:
        for i in range(0, book.nsheets):
            sheet = book.sheet_by_index(i)
            if sheet.name == sheet_name:
                break
            else:
                sheet = None
    else:
        sheet = book.sheet_by_index(sheet_index)
    if no_header:
        k = 0
        lt = list(string.uppercase)
        headers = dict((i, lt[i]) for i in range(sheet.ncols))  # todo boast lt to size of sheet.ncols
    else:
        headers = dict((i, sheet.cell_value(0, i)) for i in range(sheet.ncols))
        k = 1
    ret = list()
    for i in range(k, sheet.nrows):
        ret.append(dict((headers[j], sheet.cell_value(i, j)) for j in headers))
    return ret,headers


def xlsDictWriter(f, list_of_dicts, origin_header=None, sheet_name='sheet'):
    """"
    writes a sheet into a xls file
    :param f: string for the filename sheet_name is the name of the sheet to be written, can't be overwritten.
    :param list_of_dicts: the list of dicts to be written in a xls sheet - limited to 65536 rows due to xlwt
    :param sheet_name: string for the sheet's name to be written
    :param origin_header: dict of strings representing a (sub)set of the keys in each list_of_dicts. The desired order is designated by the order of the (numerical) keys of this dict.
    """

    #TODO: handle unicode strings - from django.utils.encoding import smart_str, smart_unicode? see sheet.write
    try:
        tmp = xlrd.open_workbook(f)
        book = copy(tmp)
    except:
        book = xlwt.Workbook(encoding="UTF-8")
    sheet = book.add_sheet(sheet_name)
    if len(list_of_dicts) > 0:
        dict_header = set()
        for row in list_of_dicts:
            for r in row.keys():
                dict_header.add(r)
        if origin_header:
            header = list()
            order = origin_header.keys()
            order.sort()
            print order
            for i in order:
                if origin_header[i] in dict_header:
                    header.append(origin_header[i])
            rest = dict_header - set(origin_header.values())
            for i in rest:
                header.append(i)
        else:
            header = dict_header
        print header
        for i, key in enumerate(header):
            sheet.write(0, i, label=key)
        for i, r in enumerate(list_of_dicts):
            for j, key in enumerate(header):
                try:
                    if type(r[key]) == int or type(r[key]) == float:
                        sheet.write(i + 1, j, label=r[key])
                    else:
                        sheet.write(i + 1, j, label=str(r[key]))
                except:
                    sheet.write(i + 1, j, label="#N/A")
        book.save(f)
        return True
    else:
        return False
