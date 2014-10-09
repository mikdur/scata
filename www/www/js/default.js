// Sätt klass 'dark' till varannan rad i tabeller
$(document).ready(function()
{
	$("table tbody tr:even").addClass('dark');
});