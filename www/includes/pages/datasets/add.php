<h1>Upload dataset</h1>

<p>Upload sequences and quality scores for the sequences as output
from the 454 instrument. Once the files are uploaded they will be
verified to ensure that the sequence names agree between the sequence
file and the quality score file.</p>

<p>If you have sequence data which is already quality-screened,
   just leave out the quality file, and that part of the quality 
   screening will be disabled. </p>

<p>Once the dataset is verified, you will recieve an email.</p>
<p><b>Uploading large datasets can take several minutes
or more. Be patient!</b></p>
<p>Fields marked with * are mandatory fields.</p>

<?php

/**
 * Verifiera först att det finns tagset redo
 * Verifiera sen att det finns parameterset att använda
 * Dessa kontroller körs individuellt så användare kan se allt som saknas
 */

$errors = Array();


// Tagsets
$tagsets = get_all_tagsets();
$num_tagsets_ready = 0;
foreach($tagsets as $tagset)
{
	$tagset = new Tagset($tagset);
	if ($tagset->is_ready())
	{
		$num_tagsets_ready++;
	}
}

/**
 * Skriv ut formuläret om inga fel inträffade
 */
if (count($errors) == 0)
{
	$form = new FormObject("");
	$form->enable_fileupload();

	$form->add_field("name", "name", "Name", "text", "mandatory=true");
	
	// Select med tagset
	$form->add_field("tagset5", "tagset5", "5' Tag set", "select", "mandatory=true");
	foreach($tagsets as $tagset)
	{
		$tagset = new Tagset($tagset);
		if ($tagset->is_owner() && $tagset->is_ready())
		{
			$form->add_option("tagset5", $tagset->id, $tagset->get_name());		
		}
	}
	$form->add_option("tagset5", "0", "none");
	$form->add_field("primer5", "primer5", "5' primer sequence. All accepted sequences are required to have the primer. The actual primer sequence is removed before clustering.", 
			 "text");

	$form->add_field("primer5score", "primer5score", "Proportional primer match (value between 0 and 1)",
			 "text", "default=0.9");

	$form->add_field("tagset3", "tagset3", "3' Tag set", "select", "mandatory=true");
	$form->add_option("tagset3", "0", "none");
	foreach($tagsets as $tagset)
	{
		$tagset = new Tagset($tagset);
		if ($tagset->is_owner() && $tagset->is_ready())
		{
			$form->add_option("tagset3", $tagset->id, $tagset->get_name());		
		}
	}
	$form->add_field("primer3", "primer3", "3' primer sequence. All accepted sequences are required to have the primer. The actual primer sequence is removed before clustering.", 
			 "text");
	$form->add_field("primer3score", "tagset3score", 
			 "Proportional primer match (value between 0 and 1)",
			 "text", "mandatory=true;default=0.9");

	$form->add_field("max_len", "max_len", 
			 "Maximum sequence length after primer and tag trimming (truncates sequences to this length starting from 5' primer, 0 implies no truncation)", "text", "default=0;mandatory=true");
	$form->add_field("min_len", "min_length",
			 "Minimum sequence length to keep (as counted after primer and tag trimming)",
			 "text", "default=200;mandatory=true");
	$form->add_field("mean_qual", "mean_qual",
			 "Minimum mean quality of bases in a read to keep",
			 "text", "mandatory=true;default=20");
	$form->add_field("min_qual", "min_qual",
			 "Minimum allowed base quality. (Drop any read containing any base with a quality below this threshold)",
			 "text", "default=10;mandatory=true");
	$form->add_field("raw_filtering", "raw_filtering", "Quality filter type. Select 'full sequence' if your data are trimmed and you want to apply the quality thresholds to the complete sequences. Select 'Extract HQR' if you want SCATA to extract the longest High Quality Region from your reads. The HQR is defined as the longest part of a read that fulfills all the quality thresholds set above. Can be good for untrimmed reads from the sequencer.", "select", "mandatory=true");
	$form->add_option("raw_filtering", "0", "Full sequences");
	$form->add_option("raw_filtering", "1", "Extract HQR");

	$form->add_field("fas_file", "fas_file",
			 "Fasta file", "file", "mandatory=true");
	$form->add_field("qual_file", "qual_file",
			 "Quality file (fasta formated), leave blank to disable filtering on quality scores", 
			 "file", "");

	// Hantera skickat formulär
	if ($form->got_post())
	{
		$form->handle_input_data();
		$dataset = new Dataset();
		$dataset->set_name($form->get_value('name'));
		$dataset->set_tagset5($form->get_value('tagset5'));
		$dataset->set_primer5($form->get_value('primer5'));
		$dataset->set_primer5score($form->get_value('primer5score'));
		$dataset->set_tagset3($form->get_value('tagset3'));
		$dataset->set_primer3($form->get_value('primer3'));
		$dataset->set_primer3score($form->get_value('primer3score'));
		$dataset->set_max_len($form->get_value('max_len'));
		$dataset->set_min_len($form->get_value('min_len'));
		$dataset->set_mean_qual($form->get_value('mean_qual'));
		$dataset->set_min_qual($form->get_value('min_qual'));
		$dataset->set_raw_filtering($form->get_value('raw_filtering'));
		$dataset->create();

		$success = true;
		// Försöker flytta uppladdade qual-filen
		if (!move_uploaded_file($_FILES["fas_file"]["tmp_name"], 
					DIR_DATASET_QUAL ."/". $dataset->id . 
					".fas"))
		  {
		    $errors[] = "Couldn't upload fasta-file.";
		    $success = false;
		  }
		// Försöker flytta uppladdade fas-filen
		move_uploaded_file($_FILES["qual_file"]["tmp_name"], 
				   DIR_DATASET_FAS ."/". $dataset->id .
				   ".qual");
		
		if ($success)
		  {
		    redirect("./?p=datasets&message=added&id={$dataset->id}");
		  }
		else
		  {
		    $dataset->delete();
		    echo "<div>";
		    echo "Error(s) occured:";
		    echo "<ul>";
		    foreach($errors as $error)
		      {
			echo "<li>{$error}</li>\n";
		      }

		  }


	}

	// Skriv ut formuläret
	echo $form->generate_javascript();
	echo $form->generate_form("./?p=datasets&amp;do=add", "Add dataset", false);
}
/**
 * Annars, skriv ut felmeddelandena
 */
else
{
	echo "Cannot add dataset:";
	echo "<ul>";
	foreach($errors as $error)
	{
		echo "<li>{$error}</li>";
	}
	echo "</ul>";
}
?>