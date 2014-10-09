<h1>Upload new reference sequence</h1>
<p>Reference sequences are used during the clustering process to
assign identities to clusters using the same criteria as
for the clustering. The reference sequences should be trimmed to
only include the sequence expected to be found among the sequenced amplicons.</p>

<p>Some sanity checks are run on the uploaded sequences. Once those checkes
are done, you will recieve an email indicating the result of the checks.</p>

<p>The fields marked with * are mandatory fields.</p>
<?php
	$success = true;
	if(isset($_POST['submit']))
	{
		$errors = Array();
		if (!isset($_POST['name']) || $_POST['name'] == "")
		{
			$errors[] = "Name not set.";
			$success = false;
		}
		if (!isset($_POST['description']) || $_POST['description'] == "")
		{
			$errors[] = "Description not set.";
			$success = false;
		}
		$cutseq = "";
		if (!isset($_POST['cutseq']) || $_POST['cutseq'] == "")
		{
		        $cutseq = "";
		}else{
		        $cutseq = $_POST['cutseq'];
		}
		$cutlen = 0;
		if (!isset($_POST['cutlen']) || $_POST['cutlen'] + 0  == 0)
		{
		        $cutlen = 0;
		}else{
		        $cutlen = $_POST['cutlen'] + 0;
		}
		if ($_FILES['file']['name'] == "")
		{
			$errors[] = "No file chosen.";
			$success = false;
		}
		if ($success == true)
		{
			$refseq = new Refseq();
			$refseq->set_name($_POST['name']);
			$refseq->set_description($_POST['description']);
			$refseq->set_cutseq($cutseq);
			$refseq->set_cutlen($cutlen);
			$refseq->set_public(isset($_POST['public']));
			$refseq->create();
			
			// Försöker flytta uppladdade filen
			if (@move_uploaded_file($_FILES["file"]["tmp_name"], DIR_REFSEQ_FAS ."/". $refseq->id . ".fas"))
			{
				$success = true;
			}
			// Gick inte att flytta, radera settet i databasen och skriv ut felmeddelande.
			else
			{
				$refseq->delete();
				$errors[] = "Couldn't upload file.";
				$success = false;
			}
		}
		if ($success)
		{
			redirect("./?p=refseqs&amp;message=added&amp;id={$refseq->id}");
		}
	}

	if (!$success)
	{
		echo "<div>";
		echo "Error(s) occured:";
		echo "<ul>";
		foreach($errors as $error)
		{
			echo "<li>{$error}</li>\n";
		}
		echo "</ul>";
		echo "</div>";
	}

?>
<form id="form_object" method="post" action="./?p=refseqs&amp;do=add" enctype="multipart/form-data">
      <div>
		<ul>
			<li>
				<label for="form_name" id="form_name_label">Name <span class="form_mandatory">*</span></label>
				<input type="text" id="form_name" name="name" />
			</li>
			<li>
				<label for="form_description" id="form_description_label">Description <span class="form_mandatory">*</span></label>
				<textarea cols="" rows="" id="form_description" name="description" ></textarea>
			</li>
<!--			<li>
				<label for="form_cutseq" id="form_cutseq_label">Cutseq - not implemented, please ignore</label>
				<input type="text" id="form_cutseq" name="cutseq" size="30"/>
			</li>
			<li>
				<label for="form_cutlen" id="form_cutlen_label">Cutlen - not implemented, please ignore</label>
				<input type="text" id="form_cutlen" name="cutlen" default="0" size="5" />
			</li>
-->
			<li>
				<label for="form_name" id="form_public_label">Check this box to make this set of reference sequences available to other users of SCATA.</label>
				<input type="checkbox" id="form_public" name="public" value="true" />
			</li>
			<li>
				<label for="form_file" id="form_file_label">File <span class="form_mandatory">*</span></label>
				<input type="file" id="form_file" value="" name="file" />
			</li>
			<li>
				<input type="submit" name="submit" id="form_submit" value="Upload file" />
			</li>
		</ul>
	</div>
</form>