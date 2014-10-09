
<?php
if (isset($_GET['id']) && $_GET['id'] != "")
{
	$parameterset = new ParameterSet($_GET['id']);
	if (!$parameterset->user_has_rights())
	{
		redirect("./?p=parametersets");
	}
	echo "<h1>Details for parameter set '" . $parameterset->get_name() . "'</h1>";
	echo "<table>";
	echo "
		<tr>
			<th>Name</th>
			<th>Description</th>
			<th>Value</th>
		</tr>";
	$parameters = get_all_parameters();
	foreach ($parameters as $parameter)
	{
		$parameter = new Parameter($parameter, $parameterset->id);
		if (!$parameter->is_available())
		{
			continue;
		}
	
		if ($parameter->get_type() == "text")
		{
			$value = $parameter->get_text_value($parameterset->id);
		}
		else
		{
			$value = $parameter->get_select_value($parameterset->id);
			$value = $value->get_value();
		}
		echo "
			<tr>
				<td><b>{$parameter->get_name()}</b></td>
				<td>{$parameter->get_description()}</td>
				<td>{$value}</td>
			</tr>\n";
	}
	echo "</table>";
}
else
{
	redirect("./?p=parametersets");
}
?>

<a href="./?p=parametersets">Back</a>