<h1>Add job</h1>
The fields marked with * are mandatory fields.
<?php

/**
 * Verifiera först att det finns dataset, tagset och referensset redo
 * Verifiera sen att det finns parameterset att använda
 * Dessa kontroller körs individuellt så användare kan se allt som saknas
 */
$errors = Array();

// Dataset
$datasets = get_all_datasets_ready();
$num_datasets_ready = 0;
foreach($datasets as $dataset)
{
	$dataset = new Dataset($dataset);
	if ($dataset->is_ready())
	{
		$num_datasets_ready++;
	}
}
if ($num_datasets_ready == 0)
{
	$errors[] = "No datasets ready.";
}

// Referenssekvenser
// Det går bra att inte skriva ut några refseqs. Ett jobb behöver inte ha något associerat refseq.
$refseqs = get_all_refseqs();
$num_refseqs_ready = 0;
foreach($refseqs as $refseq)
{
	$refseq = new Refseq($refseq);
	if ($refseq->is_ready())
	{
		$num_refseqs_ready++;
	}
}



// Parameterset
$paramsets = get_all_parametersets();
$num_parametersets_ready = count($paramsets);
if ($num_parametersets_ready == 0)
{
	$errors[] = "No parameter sets ready.";
}

/**
 * Skriv ut formuläret om inga fel inträffade
 */
if (count($errors) == 0)
{
	$form = new FormObject("");

	$form->add_field("name", "name", "Name", "text", "mandatory=true");
	//$form->add_field("description", "description", "Description", "text", "");

	// Select med parameterset
	$form->add_field("paramsets", "paramsets", "Parameter set", "select", "mandatory=true");
	foreach($paramsets as $paramset)
	{
		$paramset = new ParameterSet($paramset);
		$form->add_option("paramsets", $paramset->id, $paramset->get_name());
	}
	
	
	// Multiselect med dataset
	$form->add_field("datasets", "datasets", "Data set(s)", "select", "mandatory=true;multiselect=true");
	foreach($datasets as $dataset)
	{
		$dataset = new Dataset($dataset);
		$form->add_option("datasets", $dataset->id, $dataset->get_name());
	}

	// Multiselect med referenssekvenser
	if ($num_refseqs_ready > 0)
	{
		$form->add_field("refseqs", "refseqs", "Reference sequence(s)", "select", "multiselect=true");
		foreach($refseqs as $refseq)
		{
			$refseq = new Refseq($refseq);
			if (($refseq->is_owner() || $refseq->is_public() ) && $refseq->is_ready())
			{
				$form->add_option("refseqs", $refseq->id, $refseq->get_name());			
			}
		}
	}

	// Hantera skickat formulär
	if ($form->got_post())
	{
	  $form->handle_input_data();
	  if ($form->validate_form())
	    {

	      echo "<br/><br/>Hej" . $form->get_value('datasets') . "<br /> <br />";

	      $form->handle_input_data();
	      $job = new Job();
	      $job->set_name($form->get_value('name'));
	      //$job->set_description($form->get_value('description'));		
	      $job->set_description("Waiting for clustering to start");
	      $job->set_data_sets($form->get_value('datasets'));
	      $job->set_reference_sequences($form->get_value('refseqs'));
	      $job->set_parameter_set($form->get_value('paramsets'));
	      $job->create();
	      redirect("./?p=jobs&message=added&id={$job->id}");
	    }
	}
	// Skriv ut formuläret
	echo $form->generate_javascript();
	echo $form->generate_form("./?p=jobs&amp;do=add", "Add job", false);
}
/**
 * Annars, skriv ut felmeddelandena
 */
else
{
	echo "No jobs can be added:";
	echo "<ul>";
	foreach($errors as $error)
	{
		echo "<li>{$error}</li>";
	}
	echo "</ul>";
}
?>
