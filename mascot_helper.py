from bottle import route,get, post, run, request, template
import sqlite3 as lite
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
import os
import datetime
from os import listdir
from os.path import isfile,isdir, join
import cgi

#has to be started as process of the same user that owns mascot/sequence directory, i.e. www-data
#parse rule for the accession has to be ">[^|]*|[^|]*|\([^ ]*\)" i.e. no.34
#parse rule for the accession has to be ">[^ ]* \(.*\)" i.e. no.13


#~ import MySQLdb 
#n.b.:https://pypi.python.org/pypi/bottle-mysql !!!

# local database, used for working at home
#conn = MySQLdb.connect(host = "localhost", user = "pyramid", passwd = "pw", db = "LigandosphereDB_toy", port = 3306)
#~ conn = MySQLdb.connect(host = "192.168.123.61", user = "pyramid", passwd = "psswd", db = "LigandosphereDB_toy", port = 3306)


def init_prebuilt():
  file_sqlite = '/abi-data/walzer/immuno/mascot/prebuilt_fastas.db'
  c0 = 'CREATE TABLE fasta (name TEXT, pk INTEGER PRIMARY KEY AUTOINCREMENT, path TEXT UNIQUE, version INTEGER AUTOINCREMENT, date TEXT, url TEXT);'
  cX = list()
  names = [('swissprot_human_WOI','ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/proteomes/HUMAN.fasta.gz'),'swissprot_mouse_WOI','imgt_MAIN','imgt_ALL','imgt_c1','imgt_c2',]
  for x in names:
    #read from url - see imgt processor py
    #apply filter if neccessary
    #create dir
    #write fasta
    #get release
    #get date
    #insert to db
    cX.append('INSERT INTO fasta (sequence) VALUES (\'%s\')'% pepseq )  

#~ file_sqlite = '/tmp/mascot/prebuilt_fastas.db'
#~ c0 = 'CREATE TABLE fasta (name TEXT, pk INTEGER PRIMARY KEY AUTOINCREMENT, path TEXT UNIQUE, version INTEGER, date TEXT, url TEXT);'
#~ con = lite.connect(file_sqlite)
#~ cur = con.cursor()    
#~ cur.execute(c0)
#~ newpath = ['/tmp/mascot/prebuilt_fastas/human','/tmp/mascot/prebuilt_fastas/mouse']
#~ for i in newpath:
  #~ if not os.path.exists(i): 
    #~ os.makedirs(i)

#~ http://stackoverflow.com/a/5693214 !!!!!

#~ c1 = 'INSERT INTO fasta (name, path, date, url, version) VALUES (\'%s\',\'%s\',\'%s\',\'%s\',%s)'% ('n0',newpath[0],datetime.datetime.now(),'someurl',0)  
#~ cur.execute(c1)
#~ c2 = 'INSERT INTO fasta (name, path, date, url, version) VALUES (\'%s\',\'%s\',\'%s\',\'%s\',%s)'% ('n1',newpath[1],'now()','someurl',1 )  
#~ cur.execute(c2)

  con = lite.connect(file_sqlite)
  cur = con.cursor()    
  cur.execute(c0)
  for c in cX:
    cur.execute(c)
  con.commit()
  con.close()
  return true


@route('/')
def startpage():
  return """
    <h1>Lab&#39;s little helper</h1>
    <p></p>
    <h2>mascot:</h2>
    <ul>
      <li><a href="Mslots">slots status</a></li>
      <li><a href="Mupdate">slot update</a></li>
    </ul>
    <p>&nbsp;</p>
  """
    
@route('/hello')
def hello():
  return "Hello! But you should not be here as this site is not linked anywhere - go back to work!"
    
def get_sites():
  mypath = '/home/user/mascot/basisfasta'
  onlyfiles = [ f for f in listdir(mypath) if ( isfile(join(mypath,f)) and f.endswith('.fasta') and (not f.startswith('~')))]
  prebuilds = ['none'] + onlyfiles

  mypath = '/home/user/mascot/sequence'
  onlyslots = [ f for f in listdir(mypath) if ( isdir(join(mypath,f)) and f.startswith('slot_'))]
  slots = ['none'] + onlyslots
  
  return slots,prebuilds
  
def invalid_seq(seq):
  allowed = ['A','R','N','D','C','E','Q','G','H','I','L','K','M','F','P','S','T','W','Y','V','U','X','Z','B']
  if isinstance(seq, str) or isinstance(seq, SeqRecord):
    for s in seq:
      if s not in allowed:
        return s + ' is not allowed'
  elif isinstance(seq, list):
    for r in seq:
      for s in str(r.seq):
        if s not in allowed:
          #~ print s, 'is not allowed in', str(r.id)
          return s + ' is not allowed in ' + str(r.id)
  elif isinstance(seq, dict):
    for r in seq:
      for s in str(seq[r]):
        if s not in allowed:
          #~ print s, 'is not allowed in', str(r.id)
          return s + ' is not allowed in ' + str(r)
  else:
    return True
  #~ TODO r.id check mascot compliance  
  return False

def ambigous_aminoacid(seq):
  ret = ''
  if isinstance(seq, str) or isinstance(seq, SeqRecord):
    if 'X' in seq:
      ret = ret + 'Added sequence contains aminoacid code X!\n'
  elif isinstance(seq, list):
    for r in seq:
      if 'X' in r:
        ret = ret + 'Added ' + str(r.id) + ' contains aminoacid code X!\n'
  elif isinstance(seq, dict):
    for r in seq:
      if 'X' in seq[r]:
        ret = ret + 'Added ' + str(r) + ' contains aminoacid code X!\n'
  else:
    ret = 'Something unreadable added!'
  return ret
  
def acc_not_seen_before(what,seen):
  if what.id.split('|',1)[-1] not in seen:
    seen.add(what.id.split('|',1)[-1])
    return True
  return False  

def seq_not_seen_before(what,seen):
  if str(what.seq) not in seen:
    seen.add(str(what.seq))
    return True
  return False

def remove_dupl(seq):  
  if isinstance(seq, list):
    uni = set()
    seq = [x for x in seq if acc_not_seen_before(x,uni)]    
    uni = set()
    seq = [x for x in seq if seq_not_seen_before(x,uni)]    
  return seq

@get('/Mupdate')
def m_update():
  slots,prebuilds = get_sites()
  mupdate_form = """
    <form action="/Mupdate" method="post" enctype="multipart/form-data">
      <h2>Slotupdate:</h2>
      <table>
        <tr>
          <td>Select a slot:</td>
          <td><select name="slotselection"> %s </select></td>
        </tr><tr>
        <tr>
          <td>Select a basis:</td>
          <td><select name="basisselection"> %s </select></td>
        </tr><tr>
          <td>Choose:</td>
          <td><input type="radio" name="inputmethod" value="manual"> manual entry</td>
          <td>name:<input name="header" type="text" /> sequence:<input name="sequence" type="text" /> description:<input name="description" type="text" /></td>
        </tr><tr>
          <td>Or:</td>
          <td><input type="radio" name="inputmethod" value="file"> fasta file upload</td>
          <td><input type="file" name="fileupload" multiple /></td>
          <td>No accession:<input type="checkbox" name="uniprotonly" value="uniprotonly"> I understand the risks (and there are!)</td>
        </tr><tr>
        </tr><tr>
          <td colspan="3" align="center"><input value="Go!" type="submit", name="mupdate" style="height:40px; width:600px"/></td>
        </tr>
      </table>
    </form>
  """
  return mupdate_form % (['<option value="'+str(i)+ ('" selected="selected')[i*100:] +'">'+str(v)+'</option>' for i,v in enumerate(slots)],['<option value="'+str(i)+ ('" selected="selected')[i*100:] +'">'+str(v)+'</option>' for i,v in enumerate(prebuilds)])
    
@post('/Mupdate')
def do_m_update():
  slots,prebuilds = get_sites()
  goback = '%s Go <a href="javascript:history.back()">back.</a>'
  slot = request.forms.get('slotselection')
  if slot=="0":
    return goback % 'Select a slot first, dumdum. '
  basis = request.forms.get('basisselection')
  myfasta = list()
  baseids = list()
  records = dict()
  ### base fasta given, needs to be in swissprot header format: xx|accession|name description
  if basis!="0":  
    with open('/home/user/mascot/basisfasta/'+prebuilds[int(basis)], 'rb') as file:
      myfasta = list(SeqIO.parse(file, "fasta"))
      for r in myfasta:
        r.id= 'gln|' + r.id.split('|',1)[-1]
        r.description= 'gln|' + r.description.split('|',1)[-1]
        r.name= 'gln|' + r.name.split('|',1)[-1]
        baseids.append(r.id.split('|',1)[-1])
  met = request.forms.get('inputmethod')
  if not met:
    return goback % 'Choose a input first, buddy.'
  ent = None
  verbo = ''
  ### single entry given 
  if met == 'manual':
    ent = (request.forms.get('header').replace(" ","_"),request.forms.get('sequence').upper(),request.forms.get('description').replace('\\',"#"))
    if not ent[0] or not ent[1]:
      return goback % 'You cannot just leave something empty here, dear.'
    elif invalid_seq(ent[1]):
      return goback % 'You cannot give a non-aminoacid sequence, dear.' + invalid_seq(ent[1])
    else:
      record = SeqRecord(Seq(ent[1], IUPAC.protein), id="gnl|"+slots[int(slot)]+"|"+ent[0], name=ent[0], description=ent[2])
      verbo = ambigous_aminoacid(ent[1])
      myfasta.append(record)
  ### fasta given, needs to be in swissprot header format: xx|accession|name description
  elif met == 'file':
    uniprotonly = request.forms.get('uniprotonly') #getall for multiple files hack
    upload = request.files.getall('fileupload') #getall for multiple files hack
    if not upload:
      return goback % 'Choose a file first, buddy.'
    if not isinstance(upload,list):
      print "was no list - wtf?!"
      upload = [upload]
    i = 0
    seqset = set()
    accset = set()
    for file in upload:
      name, ext = os.path.splitext(file.filename)
      print name
      if ext not in ('.fasta', '.fa'):
        return goback % 'File extension not allowed.'
      ps = SeqIO.parse(file.file, "fasta")  
      if uniprotonly != 'uniprotonly':
        for record in ps:
            if acc_not_seen_before(record,accset) and seq_not_seen_before(record,seqset):
              records['gln|' + str(i) + 'nid|' + record.id.split('|',1)[-1]] = str(record.seq)
              i = i+1
      else:
        for record in ps:
          records[record.id+'|'+record.description] = str(record.seq)
    if uniprotonly != 'uniprotonly':
      test = invalid_seq(records)
      verbo = ambigous_aminoacid(records)
    else:
      test = False
    if test:
      return goback % 'You cannot give a non-aminoacid sequence. Meh. <br> ' + " " +test
    #~ myfasta.extend(records)
    #~ return goback % 'Not implemented yet, sorry.'
  else:
    return '<a href="javascript:history.back()">Whoopsie. (There might be a universal conspiracy going on against you.)</a>'
  try:
    mascot_inco_dir = '/home/user/mascot/sequence/' + slots[int(slot)] + "/incoming/"
    mascot_seq_dir = '/home/user/mascot/sequence/' + slots[int(slot)] + "/current/"
    version = 0
    onlyfiles = [ f for f in listdir(mascot_seq_dir) if ( isfile(join(mascot_seq_dir,f)) and f.endswith('.fasta') and (not f.startswith('~')))]
    #clean first
    mascot_old_dir = '/home/user/mascot/sequence/' + slots[int(slot)] + "/old/"
    errfiles = [ f for f in listdir(mascot_seq_dir) if ( isfile(join(mascot_seq_dir,f)) and f.endswith('.errors') and (not f.startswith('~')))]
    for errfile in errfiles:
      os.rename(mascot_seq_dir+errfile,mascot_old_dir+errfile)
      os.rename(mascot_seq_dir+errfile[:-6]+"fasta",mascot_old_dir+"defunct_"+errfile[:-6]+"fasta")
    try:
      version = max([int(x.split('.')[0].split('_')[-1]) for x in onlyfiles]) + 1
    except:  
      version  = 0
    nufile =  slots[int(slot)] + "_"+ str(datetime.datetime.now().strftime('%d%m%Y')) + "_" + str(version) + ".fasta"
    with open(mascot_inco_dir+nufile, 'w+') as open_file:
      #~ print '111 ---',mascot_seq_dir , slots[int(slot)]  , "_" , str(datetime.datetime.now().strftime('%d%m%Y'))  ,  "_"  , str(version) ,  ".fasta - length in sequences", len(myfasta) 
      SeqIO.write(myfasta, open_file, "fasta")
      for rs in records:
        r = SeqRecord(Seq(records[rs], IUPAC.protein), id=rs,description="")
        r = SeqIO.write(r, open_file, "fasta")
      os.rename(mascot_inco_dir+nufile,mascot_seq_dir+nufile)
  except:
    return '<a href="javascript:history.back()">Woah! (Cannot write into your slot.)</a>'
     
  #~ remark = ''
  #~ if id_clash:
      #~ remark = remark + 'Your fastas had conflicting ids, have been overwritten with generic ones.'
  head = ''
  try:
    with open(mascot_seq_dir+nufile, 'r') as myfile:
      N = 13*10
      head = list(islice(myfile, N))
  except Exception as e: 
    return '<a href="javascript:history.back()">Crap! (Cannot read from your slot.)</a>' + '<br> (' +  str(e)  + ')'
  
  remark = 'Your fastas had conflicting accessions, have been overwritten with generic ones.'
  return 'look for ' + nufile + ' @ mascot_status <br>' + verbo + '<br>' + remark  + '<br>'  + '<br>' + 'first N lines of the fasta on the server:' + '<br>' + '<br>'.join(head)


run(host='192.168.123.136', port=9999, reloader=True)
