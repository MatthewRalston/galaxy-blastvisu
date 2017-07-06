
<%
	import logging
	import argparse
	import sys
	import re
	from Bio.Blast import NCBIXML


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#---data necessary to the header
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	def topofmind():

		input = open(hda.file_name,"r")
		INPUT = input.read()		
		regex = re.compile('<BlastOutput_version>(.*)</BlastOutput_version>')
		version = regex.search(INPUT).group(1)
		regex = re.compile('<BlastOutput_db>(.*)</BlastOutput_db>')
		listeDb = regex.search(INPUT).group(1).split('/')
		db = listeDb[-1]

		regex = re.compile('<Iteration_query-def>(.*)</Iteration_query-def>')
		queries = regex.findall(INPUT)

		regex = re.compile('<Iteration_query-len>(.*)</Iteration_query-len>')
		length = regex.findall(INPUT)

		input.close()

		return version,db,queries,length


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#---color sequence based on the score (visualization)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	def getcolor(score):

		if score < 0 : return "white.gif"
		elif score < 40 : return "black.gif"
		elif score < 50 : return "blue.gif"
		elif score < 80 : return "green.gif"
		elif score < 200 : return "purple.gif"
		else : return "red.gif"


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#---password sequence (description)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	def sequenceName(Seq,hit_id,Db):

		if Db == "uniprot_swissprot" or Db == "uniprot" :
			ID = re.search("\|[A-Za-z0-9]*\|",hit_id).group(0).strip()
			ID = ID.replace("|","").lstrip()
			if len(Seq)>80 :
				Seq = Seq[:80]
				Seq=Seq+"..."

			return Seq,hit_id,ID,hit_id

		else :
                    # # This is overly-specific regex code that failed in every scenario i tested. It doesn't seem to do much either...
		    #     IDs = re.search("\([A-Za-z0-9\.]*\)", Seq).group(0).strip()
		    #     SeqName = Seq.replace(IDs,"",1).lstrip()
		    #     if len(SeqName)>80 :
		    #     	SeqN = SeqName[:80]
		    #     	SeqN=SeqN+"..."
		    #     else :
		    #     	SeqN = SeqName.replace(".","")	
		    #     ID = re.search("[A-Za-z0-9]+",IDs).group(0).strip()	
		    #     ote = re.search("[a-z]{2}\|[0-9]{2,}\|",hit_id).group(0).strip()
		    #     h_id = hit_id.replace( ote , ">")
		    #     hit_cmd = hit_id.replace( ote , " ")
		
		    #     return SeqN,h_id,ID,hit_cmd
                    return Seq,Seq,Seq,Seq

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#---calculation spaces for the visualization of the alignment (description)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	def espace(start1,start2):

		tai1=len(map(int,str(start1)))
		tai2=len(map(int,str(start2)))
		space1 = (10-tai1)*" "
		space2 = (10-tai2)*" "

		return space1,space2 


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#---calculation gaps for the visualization of the alignment (description)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	def gapStop(hsp,start,sbj_start,reverse,i):

		hspQ = hsp.query[i:i+60]
		gap = hspQ.count("-")
		hspS = hsp.sbjct[i:i+60]
		gap2 = hspS.count("-")
		stop1 = start+i+59-gap
		if reverse : stop2 = sbj_start-i-59+gap2
		else : stop2 = sbj_start+i+59-gap2 

		return stop1,stop2,hspQ,hspS

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#---visualization of the alignment (description)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	def align(start,end,sbj_start,sbj_stop,hsp):

		reverse = False
		TABD = ""
		gap,gap2,i = 0,0,0
		if end>start and sbj_start>sbj_stop:
			reverse= True
		space1,space2 = espace(start,sbj_start)
		stop1,stop2,hspQ,hspS = gapStop(hsp,start,sbj_start,reverse,i)
		TABD += 	"<br><br>Query %i%s%s %i" %(start,space1,hspQ,stop1)	
		TABD +=     	"<br>                %s" %hsp.match[i:i+60]
		TABD +=     	"<br>Sbjct %i%s%s %i" %(sbj_start,space2,hspS,stop2)	
		i += 60
		while i+start < end :
			start1=stop1+1
			if reverse : start2=stop2-1 
			else :	start2=stop2+1 
			space1,space2 = espace(start1,start2)
			stop1,stop2,hspQ,hspS = gapStop(hsp,start,sbj_start,reverse,i)
			if stop1>end : stop1 = end
			if reverse :														
				if stop2<sbj_stop : stop2 = sbj_stop
			else:
				if stop2>sbj_stop : stop2 = sbj_stop 
			TABD += "<br><br>Query %i%s%s %i " %(start1,space1,hspQ,stop1)
			TABD +=     "<br>                %s" %hsp.match[i:i+60]
			TABD +=     "<br>Sbjct %i%s%s %i" %(start2,space2,hspS,stop2)
			i += 60
		TABD += "<br><br>"

		return TABD

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#---description of the alignment (description)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	def descri(hsp,Seq,hit_id,Db,Name,end,start):
	
		TABD=""
		iden = hsp.identities
		score	= hsp.bits
		evalue = hsp.expect
		sbj_start = hsp.sbjct_start
		sbj_stop = hsp.sbjct_end					
		SeqN,h_id,ID,hit_cmd = sequenceName(Seq,hit_id,Db)
		TABD += "<tr><td><input type='checkbox' onclick='change(this.id)' name='list' id='BoxQuery%sSeq%s' value=%i seq='%s\n'>" % (Name,h_id,score,hit_cmd)
		TABD += "<a href='javascript: toggleDisplay(\"Query%sSeq%s\");' id ='Query%s%s' >%s</a>" % (Name,h_id,Name,h_id,SeqN)
		TABD += "</input></td>"
		TABD += "<td align='center'>%s</td>" %score
		TABD +=	"<td align='center'>%s</td>" %evalue
		if Db in SITES.keys() :
			TABD +=	"<td align='center' ><a href='%s/%s' target='_blank'><u>%s</u></a></td>" % (SITES[Db],ID,ID)
		else : 
			TABD += "<td align='center'>%s</td>" %ID
		TABD += "</tr><tr>"
		TABD += "<td style='display: none; margin-left:20px; border:1px solid #999999;' id='Query%sSeq%s' colsan='4'>" %(Name,h_id)
		TABD += "<pre style='font-size:12px;line-height:1.42857;'>" 
		TABD += "<br>%s%s<br>Length=%i" %(h_id,Seq,alignment.length)
		TABD +=	"<br> Score = %i bits , Expect = %i" %(score,evalue)
		TABD += "<br> Identities = %i/%i , Gaps = %i/%i" %(iden,hsp.align_length,hsp.gaps,hsp.align_length)
		TABD += align(start,end,sbj_start,sbj_stop,hsp)
			
		return TABD

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#---plotting of multiple hsp (visualization)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	def hspgraph(indice,iteration,begin,start,stop,end,strip,col,lg,width,color):

		TABV=""
		while indice <= iteration :
			if begin[indice]<start :
				if begin[indice] != 1 :	TABV += "<img width='%i' height='2' src='/static/visuBlast/Images/white.gif'/>" % (500*begin[indice]/lg)	
				if stop[indice]<start or stop[indice] == start:
					TABV += "<img width='%i' height='4' src='/static/visuBlast/Images/%s'/>" % (strip[indice],col[indice])
					if stop[indice] == start :	TABV += "|"
					else : TABV += "<img width='%i' height='2' src='/static/visuBlast/Images/white.gif'/>" % (500*(start-stop[indice])/lg)
				else :	
					TABV += "<img width='%i' height='4' src='/static/visuBlast/Images/%s'/>" % (500*(start-begin[indice]+1)/lg,col[indice])
					TABV += "|"
				TABV += "<img width='%i' height='4' src='/static/visuBlast/Images/%s'/>" % (width,color)	
			else :
				if start != 1 :	TABV += "<img width='%i' height='2' src='/static/visuBlast/Images/white.gif'/>" % (500*start/lg)		
				TABV += "<img width='%i' height='4' src='/static/visuBlast/Images/%s'/>" % (width,color)
				if end == begin[indice] :	
					TABV += "|"
					TABV += "<img width='%i' height='4' src='/static/visuBlast/Images/%s'/>" % (strip[indice],col[indice])	
				elif end > begin[indice] : 
					TABV += "|"
					TABV += "<img width='%i' height='4' src='/static/visuBlast/Images/%s'/>" % (500*(stop[indice]-end+1)/lg,col[indice])
				else : 	
					TABV += "<img width='%i' height='2' src='/static/visuBlast/Images/white.gif'/>" % (500*(begin[indice]-end)/lg)
					TABV += "<img width='%i' height='4' src='/static/visuBlast/Images/%s'/>" % (strip[indice],col[indice])	
			indice += 1	

		return TABV

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#---
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	Version,Db,Queries,Length = topofmind()

	SITES = {"uniprot" : "http://www.uniprot.org/uniprot",
	"uniprot_swissprot" : "http://www.uniprot.org/uniprot",
	"nr" : "http://www.ncbi.nlm.nih.gov/protein",
	"genbank" : "http://www.ncbi.nlm.nih.gov/nucleotide",
	"nt" : "http://www.ncbi.nlm.nih.gov/nucleotide"}
	result_handle = open(hda.file_name)
	blast_records = NCBIXML.parse(result_handle)
	Name = ""
	compt = 0
	TABV,TABD = {},{}
	for blast_record in blast_records:
		lg = blast_record.query_letters
		for alignment in blast_record.alignments:
			hit_id = alignment.hit_id
			Seq = alignment.hit_def
			nbhsp = 0
			hspmu = False
			for hsp in alignment.hsps:
				nbhsp += 1
			if nbhsp != 1 : hspmu = True
			for hsp in alignment.hsps:
				start = hsp.query_start
				end = hsp.query_end

#tab Visu & Descri
				if Queries.index(blast_record.query) != Name : 
					compt = 0
					Name = Queries.index(blast_record.query)
					TABV[Name] = "<table border='0' cellpadding='0' cellspacing='0'>"
					TABD[Name] = ""
#tabV
				if compt == 50 : 
					TABV[Name] += "</table>"
				elif compt < 50 :
					breadth = end-start+1
					color = getcolor(hsp.bits)
					width = 500*breadth/lg
					TABV[Name] += "<tr><td/><img height='2' src='/static/visuBlast/Images/white.gif' width='50'/>"
					TABV[Name] += "<tr><td/><img height='4' src='/static/visuBlast/Images/white.gif' width='50'/><td>"
					TABV[Name] += "<a href = '#Query%s' id = 'linkancre' title = '%s'>" % (hit_id,hit_id)
#tabD
				if compt == 250 :	
					TABD[Name] += "</table>"
				elif compt < 250 :
					compt += 1
					TABD[Name] += descri(hsp,Seq,hit_id,Db,Name,end,start)

#if several alignment with same sequence
					if hspmu :
						iteration,indice = 1,2
						begin,stop,strip,col = {},{},{},{}
						for hsp in alignment.hsps:
							if hsp.query_start == start : continue			
							else :
								iteration += 1
								start2 = hsp.query_start
								end2 = hsp.query_end
								begin[iteration] = hsp.query_start
								stop[iteration] = hsp.query_end
								strip[iteration] = 500*(stop[iteration]-begin[iteration]+1)/lg
								col[iteration] = getcolor(hsp.bits)
#tabV
						TABV[Name] += hspgraph(indice,iteration,begin,start,stop,end,strip,col,lg,width,color)	
#tabD
						TABD[Name] +=	"<br> Score = %i bits , Expect = %i" %(hsp.bits,hsp.expect)
						TABD[Name] += "<br> Identities = %i/%i , Gaps = %i/%i" %(hsp.identities,hsp.align_length,hsp.gaps,hsp.align_length)
						TABD[Name] += align(start2,end2,hsp.sbjct_start,hsp.sbjct_end,hsp)
						TABD[Name] += "</pre><br/></td></tr>"
#tabV
					else :
						if compt < 50 :
							if start != 1 :	TABV[Name] += "<img width='%i' height='2' src='/static/visuBlast/Images/white.gif'/>" % (500*start/lg)		
							TABV[Name] += "<img width='%i' height='4' src='/static/visuBlast/Images/%s'/>" % (width,color)
					TABV[Name] += "</a></td></tr>"					
#tabD
					TABD[Name] += "</pre><br/></td></tr>"
					break

%>




##-----------------------------------------------------------------------------------------#
##------------------------------TEMPLATE HTML CREATE---------------------------------------#


<!DOCTYPE HTML>
<html>
	<head>
		<meta http-equiv="Content-Type" content="text/html; charset=utf-8" />
   	<title>Visualization Output Blast</title>
		<link type="text/css" rel="Stylesheet" media="screen" href="/static/style/base.css">
		<link type="text/css" rel="Stylesheet" media="screen" href="/plugins/visualizations/blastvisu/static/style.css">
		<script type="text/javascript" src="/plugins/visualizations/blastvisu/static/visublast.js"></script>
	</head>
	<body style = "min-width: 900px;">
		<div class="toolForm">
		<%
				TITLE = "ABIMS : %s<br/>Database : %s<br/>" % (Version,Db)
				if len(Queries) > 1 : 
					TITLE += "Query : <select id='selectQuery' onchange='ChangeQuery()' >"
					for q,query in enumerate(Queries) : TITLE += "<option value='%s'>%s</option>" % (q,query)
					TITLE += "</select><br/>"
				else :	TITLE += "Query : %s<br/>" % Queries[0]
				TITLE += "Length : <span id = 'length'>%s</span>" % Length[0]
			%>
			<div class="toolFormTitle">${TITLE}</div>

		% for q,query in enumerate(Queries) :
				<%
					L = Length[q]
				%>
				% if  q == 0 :
					<div id="Contener_Query${q}"  length="${L}">
			 	% else :
					<div id="Contener_Query${q}"  length="${L}" style="display:none;">
				%endif

##~~~~~~~SEQUENCES AREA~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			<form>
				<textarea id="DownloadZone_Query${q}" readonly="readonly" wrap="off" class="textarea" style="display: none;">""</textarea>
			</form>	

##~~~~~~~FORM OF SEQUENCES~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			<form name="myform" id="FormList_Query${q}">

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##~~~~~~~~~~~~~~~~~~~~~VISUALISATION~~~~~~~~~~~~~~~~~~~~~~

				<div class="form-row">	
					<h3 class="form-row form-actions">		
						<a href="javascript: reduce('visu_query${q}','cadre_query${q}');">
							<span><i class="fa fa-chevron-circle-up" id="visu_query${q}"></i>	Visualization : </span>
						</a>	
					</h3>
				</div>
				<div id="cadre_query${q}" class="cadre">
					<table id="visualization_query${q}" class="visualization">
 						<tr><td align=right><img width="500" height="40" src="/static/visuBlast/Images/score.gif"/></td></tr>
						<tr><td/><img width="550" height="4" src="/static/visuBlast/Images/white.gif"/></tr>
						<tr><td/><img width="550" height="10" src="/static/visuBlast/Images/query_no_scale.gif"/></tr>
					  <tr><td/><img width="550" height="4" src="/static/visuBlast/Images/white.gif"/></tr>
						<tr><td/><img width="550" height="4" src="/static/visuBlast/Images/white.gif"/></tr>
						${TABV[q]}
					</table>
				</div>

##~~~~~~~~~~~~~~~~~~~~~~~~~TOOLS~~~~~~~~~~~~~~~~~~~~~~~~
				<div class="form-row">
					<h3 class="form-row form-actions">
						<a href="javascript: reduce('tool_query${q}','divbutton_query${q}');">
							<span><i class="fa fa-chevron-circle-up" id="tool_query${q}"></i>	Tools : </span>
						</a>						
					</h3>
					<div id="divbutton_query${q}" class="form-row divbutton">
						<input type="text" value="160" id="ScoreField_query${q}" size="7">
						<input type="button" value="Score" onClick="score(this.form.list, document.getElementById('ScoreField_query${q}').value)"> 
						<br>
						<input type="text" value="10" id="BestHits_query${q}" size="7">
						<input type="button" value="CheckBest" onClick="besthits(this.form.list, document.getElementById('BestHits_query${q}').value)">
						<br>
						<input type="text" value="Things" id="KeyWords_query${q}" size="7">
						<input type="button" value="FindWord" onClick="keywords(this.form.list,document.getElementById('KeyWords_query${q}').value)"> 
						<br>	
					</div>
				</div>	

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##~~~~~~~~~~~~~~~~~~~~~~~~~DESCRIPTIONS~~~~~~~~~~~~~~~~~~~~~~~~

				<div class="form-row">
					<h3 class="form-row form-actions">
						<a href="javascript: reduce('descript_query${q}','table_query${q}');">
							<span><i class="fa fa-chevron-circle-up" id="descript_query${q}"></i>	Descriptions : </span>
						</a>
					</h3>
				</div>
				<div id="table_query${q}">
					<table cellpadding="2" rules="all" cellspacing="0" class="description">
						<col width="70%"/>
						<col width="15%"/>
						<col width="15%"/>
						<col width="20%"/>	
						<tr><th>

##~~~~~~~~~~~~~~~~~~~~~~~~~CASEMENU~~~~~~~~~~~~~~~~~~~~~~~
							<div id="casemenu_Query${q}" class="casemenu">
								<a class="casemenubut">
									<span onClick="check(document.getElementById('FormList_Query${q}').list)" class="fa fa-check-square-o" title="Check all"> 
									</span>
								</a>
								<a class="casemenubut">
									<span onClick="uncheck(document.getElementById('FormList_Query${q}').list)" class="fa fa-external-link" title="Uncheck all"></span>
								</a>
								<a class="casemenubut aButDownOn" id="aButDownOn_Query${q}" style="border-right:1px solid #BFBFBF;">
									<span class="fa fa-floppy-o" title="Download">
									</span>
								</a>
							</div>
							<div style="margin-top:7px;">
								<p id="NumSelected_Query${q}" class="form-row NumSelected">Selected : 0</p>
							</div>
							<div style="margin-top:10px;margin-bottom:0px;margin-left:5px;">
								<text>Descriptions</text>
							</div>
							<th><text>Score</text></th>
							<th><text>Evalues</text></th>
							<th><text>Accessions</text></th>
							${TABD[q]}																										
					</div>																														
				</form>																																
			</div>																																	
			% endfor
		</div> 																																		
	</body>
</html>

