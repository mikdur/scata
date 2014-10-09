<h1>Upload Tagset</h1>
<p>The tagset defines what tags have been used to identify samples and provides a link
between the tag sequence and a symbolic name for the tag set.</p>
<p>The format is one tag per line, with tag name and tag sequence separated by a semicolon.
The length of the tag sequence has to be the same for all tags. A tag set might look like
this:
<pre>Tag1;AACTGG
Tag2;ACGATG
Tag3;TATGAA
</pre>
</p>
<p>Once the tag set is uploaded it will be validated and you will recieve an email once this is done.</p>

<p>The fields marked with * are mandatory fields.</p>

<?php
	if(isset($_POST['submit']))
	{
		print_r($_POST['submit']);
		if (!isset($_POST['name']) || !isset($_FILES['file']))
		{
			$message = "Error: Name and/or file not set";
			$success = false;
		}
		else
		{
			$tagset = new Tagset();
			$tagset->set_name($_POST['name']);
			$tagset->create();
			
			// Försöker flytta uppladdade filen
			if (@move_uploaded_file($_FILES["file"]["tmp_name"], DIR_TAGSET_TXT ."/". $tagset->id . ".txt"))
			{
				$success = true;
			}
			// Gick inte att flytta, radera settet i databasen och skriv ut felmeddelande.
			else
			{
				$tagset->delete();
				$message = "Couldn't move file.";
				$success = false;
				$tagset->delete();
			}
		}
		if ($success)
		{
			redirect("./?p=tagsets&message=added&id={$tagset->id}");
		}
		else
		{
			echo "<div>{$message}</div>";
		}
	}
?>
<form id="form_object" method="post" action="./?p=tagsets&amp;do=add" enctype="multipart/form-data">
	<div>
		<ul>
			<li>
				<label for="form_name" id="form_name_label">Name <span class="form_mandatory">*</span></label>
				<input type="text" id="form_name" name="name" />
			</li>
			<li>
				<label for="form_file" id="form_file_label">File <span class="form_mandatory">*</span></label>
				<input type="file" id="form_file" value="" name="file" />
			</li>
			<li>
				<input type="submit" name="submit" id="form_submit" value="add" />
			</li>
		</ul>
	</div>
</form>