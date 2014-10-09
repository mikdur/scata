<h1>Define new parameter set</h1>

<?php
$form = new FormObject("");
$parameter = new Parameter();
$user = new User();

$form->add_group("general", "General");
$form->add_group("parameter", "Parameter specifics");
$form->add_field("name", "name", "<b>Name</b>", "text", "mandatory=true", "general");
$form->add_field("description", "description", "<b>Description</b>", "textarea", "mandatory=true", "general");

$parameters = get_all_parameters();
for ($i = 0; $i < count($parameters); $i++)
{
	$parameter = new Parameter($parameters[$i]);

	if (!$parameter->is_available())
	{
		continue;
	}

	// Skapa text för labeln
	$label_text = "
	<b>{$parameter->get_name()}</b><br />
	{$parameter->get_description()}";
	
	if ($parameter->get_type() == "text")
	{
		$form->add_field("param-".$parameter->id, "param-{$parameter->id}", $label_text, "text", "mandatory=true", "parameter");
		$form->set_value("param-".$parameter->id, $parameter->get_default_value()->get_value());
	}
	elseif ($parameter->get_type() == "select")
	{
		$form->add_field("param-".$parameter->id, "param-{$parameter->id}", $label_text, "select", "mandatory=true", "parameter");
		$parametervalues = $parameter->get_possible_select_values();
		foreach ($parametervalues as $value)
		{
			$form->add_option("param-".$parameter->id, $value->id, $value->get_display_value(), $value->is_default());
		}
	}
	else
	{}
}

// Hantera skickat formulär
if ($form->got_post())
{
	$form->handle_input_data();
	if ($form->validate_form())
	{
		$parameterset = new ParameterSet();
		$parameterset->set_name($form->get_value('name'));
		$parameterset->set_description($form->get_value('description'));
		$parameterset->set_owner($user->get_id());
		$parameterset->write_to_db();
		
		foreach ($parameters as $parameter)
		{
			$parameter = new Parameter($parameter, $parameterset->id);
			if ($parameter->get_type() == "text")
			{
				$parameter->set_text_value($form->get_value("param-{$parameter->id}"));
			}
			elseif ($parameter->get_type() == "select")
			{
				$parameter->set_select_value($form->get_value("param-{$parameter->id}"));
			}
			else
			{
				throw new Exception("Add: Incorrect parametertype");
			}
			
		}
		redirect("./?p=parametersets&message=added&id={$parameterset->id}");
	}
}

echo $form->generate_javascript();
echo $form->generate_form("./?p=parametersets&amp;do=add", "Add set", false);

?>
