<?php
if (isset($_GET['message']))
{
	$messages = Array();
	$messages['deleted'] = "Data set has been deleted";
	$messages['added'] = "Added data set";
	notify($messages[$_GET['message']]);
}
?>

<h1>Datasets</h1>
<span class="medium"><a href="./?p=datasets&amp;do=add">Upload dataset</a></span><br />
<?php

echo "
<table>
	<thead>
		<tr>
			<th>Dataset Name</th><th>Description</th><th>Upload date</th><th>Status</th><th>Delete</th>
		</tr>
	</thead>
	<tbody>";
	
$datasets = get_all_datasets();
$form_num = 0;
if (count($datasets) > 0)
{
	foreach($datasets as $dataset)
	{
		$dataset = new dataset($dataset);
		echo "
			<tr>
				<td valign='top'>{$dataset->get_name()}</td>
                                <td><pre>{$dataset->get_description()}</pre></td>
				<td valign='top'>".$dataset->get_created_date()."</td>
				<td>".($dataset->is_ready() ? "<b>OK</b>" : "pending validation")."</td>
				<td valign='top'>
					<form method=\"post\" action=\"./?p=datasets&do=delete\" id=\"form_{$form_num}\">
						<input type=\"hidden\" value=\"{$dataset->id}\" name=\"id\">
						<a href=\"#\" onClick=\"confirm_delete({$form_num});\"><img src=\"./images/delete.png\" alt=\"delete\" /></a>
					</form>
				</td>
			</tr>";
			$form_num++;
	}
}
else
{
	echo "<tr><td>No datasets added</td></tr>";
}
echo "
	</tbody>
</table>";
?>
