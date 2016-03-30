	//--------------------
	// Globales Variables
	//--------------------

	var link = "false";
	var pinit = "false";
	var query = "0";

	//-----------------------------------------
	// Complete textarea to download sequences 
	//-----------------------------------------

	function addfile(seq)
	{ 
		var text = document.getElementById("DownloadZone_Query" + query);
		text.value = text.value + seq;
	}

	function updateSequences()
	{
		var text = document.getElementById("DownloadZone_Query" + query);
		text.value = "";
		var field = document.getElementById('FormList_Query' + query).list;
		for (i = 0; i < field.length; i++) 
		{
			if (field[i].checked == true) {addfile(field[i].getAttribute("seq"));}	
		}
	}

	//-------------------------------------
	// Control sequence/alignement display
	//-------------------------------------

	function toggleDisplay(elmt)
	{
		elmt = document.getElementById(elmt);
	  if(elmt.style.display == "none")
		{
		  elmt.style.display = "";
	  } else { 
		  elmt.style.display = "none";
		}
	}

	function reduce(id,elmt)
	{
		id = document.getElementById(id);
		elmt = document.getElementById(elmt);
		if(elmt.style.display == "none")
		{
			elmt.style.display = "";
			id.className = "fa fa-chevron-circle-up";
		} else {
			elmt.style.display = "none";
			id.className = "fa fa-chevron-circle-down";
		}
	}

	//-----------------------------------------
	// Update the number of selected sequences
	//-----------------------------------------

	function updateSelected(num)
	{
		var p = document.getElementById("NumSelected_Query" + query);
		p.innerHTML = "Selected : " + num;
		p.value = num;
		if ((num > 0) && (link =="false")) { download("true"); }
	}

	function change(id)
	{	
		var n = getselected();
		var box = document.getElementById(id);
		if (box.checked == true)
		{
			addfile(box.getAttribute("seq"));
			updateSelected(n+1);
			download("true");	
		} else {
			updateSequences();
			updateSelected(n-1);
			if (n == 1) { download("false"); }
			if (n > 1) { download("true"); }
		}		
	}

	function getselected()
	{
		var p = document.getElementById("NumSelected_Query" + query);
		if (pinit == "false") 
		{ 
			p.value = 0;
			pinit = "true";
		}
		var n = parseInt(p.value,10);
		return n;
	}

	//-----------------------------------
	// Control display of download link
	//-----------------------------------

	function download(linkwant)
	{
		if (linkwant == "true") 
		{
			//document.getElementById("ButDownOff").setAttribute("id", "ButDownOn");
			//var a = document.getElementById("divbutton").appendChild(document.createElement("a")); 
			var a = document.getElementById("aButDownOn_Query" + query);
			a.download = "idseq.txt";
			a.href = "data:text/txt," + escape(document.getElementById("DownloadZone_Query" + query).value); //???????????????????????null????????????????????
			link ="true";
		}

		if (linkwant == "false") 
		{
			//document.getElementById("ButDownOn").setAttribute("id", "ButDownOff");
			var a = document.getElementById("aButDownOn_Query" + query);
			a.removeAttribute("href") 
			link ="false";
			var text = document.getElementById("DownloadZone_Query" + query);
			text.value = "";
		}
	}

	//-------------------------------
	// Permit to check all the boxes
	//-------------------------------

	function check(field) 
	{
		var text = document.getElementById("DownloadZone_Query" + query);
		text.value = "";

		if (typeof field.length == 'undefined')
		{
			field.checked = true;
			addfile(field.getAttribute("seq"));
			updateSelected(1);
		} else {
			for (i = 0; i < field.length; i++) 
			{
				field[i].checked = true;
				addfile(field[i].getAttribute("seq"));	
			}	
			updateSelected(field.length);
		}

		download("true");
	}

	//---------------------------------
	// Permit to uncheck all the boxes
	//---------------------------------

	function uncheck(field) 
	{
		var text = document.getElementById("DownloadZone_Query" + query);
		text.value = "";
		updateSelected(0);
		if (link == "true") { download("false"); }

		if (typeof field.length == 'undefined')
		{
			field.checked = false;
			return
		}

		for (i = 0; i < field.length; i++)
		{
			field[i].checked = false;
			
		}
	}

	//-----------------------------------
	// Permit to check the n first boxes
	//-----------------------------------

	function besthits(field,num) 
	{
		if (typeof field.length == 'undefined')
		{
			if (num >= 1) { field.checked = true; }
			return
		}
		uncheck(field);
		if (num < field.length) 
		{ 
			var toNum = num;
		} else if (num >= field.length) {
			var toNum = field.length; 
		} else {
			return False;
		}
		for (i = 0; i < toNum; i++) 
		{
			field[i].checked = true;
			addfile(field[i].getAttribute("seq"));
		}
		updateSelected(toNum);
		download("true");
		
	}
	
	//--------------------------------------------------
	// Permit to check all the boxes with better score
	//--------------------------------------------------

	function score(field,value)
	{
		if (typeof field.length == 'undefined')
		{
			if (parseInt(field.value,10) >= parseInt(value,10)) { field.checked = true; }
			return
		}
		uncheck(field);
		var count = 0;
		for (i = 0; i < field.length; i++)
		{
			if (parseInt(field[i].value,10) >= parseInt(value,10))
			{
				field[i].checked = true;
				addfile(field[i].getAttribute("seq"));
				count++;
			}
		}
		updateSelected(count);
		download("true");
	}

	//---------------------------------------------------------
	// Permit to check all the boxes corresponding to keywords
	//---------------------------------------------------------

	function keywords(field,words)
	{
		if (typeof field.length == 'undefined')
		{
			if (field.getAttribute("seq").toLowerCase().indexOf(words.toLowerCase()) > -1) { field.checked = true; }
			return
		}
		uncheck(field); 
		var count = 0;
		for (i = 0; i < field.length; i++)
		{
			if ( field[i].getAttribute("seq").toLowerCase().indexOf(words.toLowerCase()) > -1)
			{
				field[i].checked = true;
				addfile(field[i].getAttribute("seq"));
				count++;
			}
		}
		updateSelected(count);
		download("true");
	}

	//---------------------------------------------
	// Permit to check up to a field already check
	//---------------------------------------------

	function upto(field)
	{
		var n = getselected();
		var count = 0;
		for (i = 0; i < field.length; i++)
		{
			if (field[i].checked == true) { break ;}
			field[i].checked = true;
			addfile(field[i].getAttribute("seq"));
			count++;
		}
		updateSelected(count + n);
		download("true");
	}


	//----------------------------------------------
	// Change display in case of multi fasta request
	//----------------------------------------------

	function ChangeQuery()
	{
		var listBox = document.getElementById("selectQuery");
		var v = listBox.options[listBox.selectedIndex].value;
		var l = document.getElementById("length");

		l.innerHTML = document.getElementById("Contener_Query"+v).getAttribute("length");   
		pinit = "false";

		toggleDisplay("Contener_Query"+query);
		toggleDisplay("Contener_Query"+v);

		query = v;
	}
