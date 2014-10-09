<?php
if (isset($_GET['message']))
{
	$messages = Array();
	$messages['deleted'] = "Parameter set has been deleted";
	$messages['added'] = "Added parameter set";
	notify($messages[$_GET['message']]);
}
?>

<h1>Parameter Sets</h1>
<span class="medium"><a href="./?p=parametersets&amp;do=add">Add parameter set</a></span><br />
<?php
if (isset($_GET['message']))
{
	$messages = Array();
	$messages["added"] = "Parameter set added.";
	
	if (array_key_exists($_GET['message'], $messages))
	{
		echo $messages[$_GET['message']];
	}
}

echo "
<table>
	<thead>
		<tr>
			<th>Name</th><th>Description</th>
		</tr>
	</thead>
	<tbody>";

$parametersets = get_all_parametersets();
$form_num = 0;
if (count($parametersets) > 0)
{
	foreach ($parametersets as $parameterset)
	{
	$parameterset = new ParameterSet($parameterset);
	echo "
		<tr id=\"paramset-{$parameterset->id}\">
			<td><a href=\"./?p=parametersets&amp;do=details&amp;id={$parameterset->id}\">".$parameterset->get_name()."</a></td>
			<td>".$parameterset->get_description()."</td>
			<td>
				<form method=\"post\" action=\"./?p=parametersets&do=delete\" id=\"form_{$form_num}\">
					<input type=\"hidden\" value=\"{$parameterset->id}\" name=\"id\">
					<a href=\"#\" onClick=\"confirm_delete({$form_num});\"><img src=\"./images/delete.png\" alt=\"delete\" /></a>
				</form>
			</td>
		</tr>";
		$form_num++;
	}
}
else
{
	echo "<tr><td>No parameter sets added</td></tr>";
}
echo "
	</tbody>
</table>";
?>
