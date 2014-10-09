<h1>Login to SCATA</h1>

<?php
if (isset($_POST['login']))
{
	$user = new User();
	$login_check = $user->login($_POST['email'], $_POST['password']);
	if ($login_check == USER_SUCCESS)
	{
		redirect("./?p=home");
	}
	else
	{
		echo "<h2>Error</h2>";
		switch ($login_check)
		{
			case USER_INVALID_CREDENTIALS :
				$error_message = "Invalid user credentials.";
				break;
			case USER_NOT_ENABLED :
				$error_message = "User is not enabled.";
				break;
			default :
				$error_message = "";
				break;
		}
		echo "<p>{$error_message}</p>";
	}
}
?>

<form method="post" action="./?p=login" id="form_object">
	<div>
		<ul>
<? if($_SERVER['HTTPS'] != "on"){
   echo "If you prefer to use a secure connection to Scata, follow this link: ".
   "<a href='https://scata.mykopat.slu.se'>".
   "https://scata.mykopat.slu.se/</a>";
 }
?>
			<li>
				<label for="email">
					E-mail address
   <? if($_SERVER['HTTPS'] == "on"){
					  echo "or SLU AD username.";
					}
?>
				</label>
				<input type="text" name="email" id="email" />
			</li>
			<li>
				<label for="password">
					Password
				</label>
				<input type="password" name="password" id="password" />
			</li>
			<li>
				<input type="submit" value="Login" name="login" id="login" />
			</li>
		 <li>
		<a href="./?p=lost_password">Lost password?</a>
</li>
<? if($_SERVER['HTTPS'] == "on"){
   echo "<li>You are using a secure connection to Scata. If you are a SLU ".
   "user, please consider using your AD username and password to log in ".
   "without the need to register. Just type your username without AD\\ in ".
   "front of it.</li>";
 }else{
   echo "<li>Please consider using ".
     "<a href='https://scata.mykopat.slu.se/'>".
     "https://scata.mykopat.slu.se/</a> instead. ".
     "That will give you a secure connection to Scata, where all informations ".
     "will be encrypted before being sent over the net. ".
     "For SLU users this will provide ".
     "the possibility to log in using your AD username and password as well. </li>";
 }
?>
</ul>
	</div>
</form>