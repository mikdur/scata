<?php
if (isset($_GET['message']))
{
	$messages = Array();
	$messages['deleted'] = "Tag set has been deleted";
	$messages['added'] = "Added tag set";
	notify($messages[$_GET['message']]);
}
?>

<h1>Tagsets</h1>
<span class="medium"><a href="./?p=tagsets&amp;do=add">Upload tag set</a></span><br />
<?php

echo "
<table>
	<thead>
		<tr>
			<th>Name</th><th>Status</th><th>Delete</th>
		</tr>
	</thead>
	<tbody>";
	
$tagsets = get_all_tagsets();
$form_num = 0;
if (count($tagsets) > 0)
{
	foreach($tagsets as $tagset)
	{
		$tagset = new Tagset($tagset);
		echo "
			<tr>
				<td>{$tagset->get_name()}</td>
				<td>".($tagset->is_ready() ? "<b>OK</b>" : "pending validation")."</td>
				<td>
					<form method=\"post\" action=\"./?p=tagsets&do=delete\" id=\"form_{$form_num}\">
						<input type=\"hidden\" value=\"{$tagset->id}\" name=\"id\">
						<a href=\"#\" onClick=\"confirm_delete({$form_num});\"><img src=\"./images/delete.png\" alt=\"delete\" /></a>
					</form>
				</td>
			</tr>";
		$form_num++;
	}
}
else
{
	echo "<tr><td>No tag sets added</td></tr>";
}
echo "
	</tbody>
</table>";
?>
