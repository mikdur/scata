<h1>Register new account</h1>
<?php
// Aktivering av konto
if (isset($_GET['do']) && $_GET['do'] == "activate")
{
		echo "<h2>Account activation</h2>";
	if (!isset($_GET['email']) && !isset($_GET['key']))
	{
		echo "
		<p>
			E-mail and key must be specified. Please use the entire link to activate.
		</p>";
		exit(0);
	}
	else
	{
		$registration = new Registration(0,$_GET['email']);

		if ($registration->activate_user($_GET['key']))
		{
			echo "
			<p>
				Your account has successfully been activated. Proceed to <a href=\"./?p=login\">login page</a>.
			</p>";
			exit(0);
		}
		else
		{
			echo "
			<p>
                An error occurred when activating the account.<br />
			</p>";
			exit(0);
		}
	}
	
}

// Hantera skickat formulär
if (isset($_POST['register']))
{
	$registration = new Registration();
	$email = Array($_POST['email_1'], $_POST['email_2']);
	$password = Array($_POST['pass_1'], $_POST['pass_2']);
	$secret = $_POST['foo'];
	$error_message = Array();
	
	
	// Se om e-post matchar och ej redan finns
	$email_check = $registration->check_email($email[0], $email[1]);
	switch ($email_check) {
		case REGISTRATION_EMAIL_EXISTS :
			$error_message[] = "User already exists.";
			break;
		case REGISTRATION_NOT_MATCHING_EMAILS :
			$error_message[] = "E-mail addresses don't match.";
			break;
		default :
			break;
	}
	
	// Se om lösenord är korrekt angivet
	$password_check = $registration->check_passwords($password[0], $password[1]);
	switch ($password_check)
	{
		case REGISTRATION_NOT_MATCHING_PASSWORDS :
			$error_message[] = "Passwords don't match.";
			break;
		case REGISTRATION_NOT_VALID_PASSWORD :
			$error_message[] = "Password is invalid.";
			break;
		default :
			break;
	}
	if (count($error_message) == 0 && $secret == "")
	{
		try
		{
			$registration->register_user($email, $password);
			$message = "Your account has been successfully registered. Please check your e-mail inbox for the activation mail.";
		}
		catch (Exception $e)
		{
			echo "Error occurred: " . $e->getMessage() . "\n";
		}		
	}
	else
	{
	}
}	
?>
<?php
if (isset($message))
{
	echo "<h2>Message</h2>";
	// Finns det något felmeddelande?
	if (count($error_message) > 0)
	{
		echo "<p>
		Error:
			<ul>";
		foreach($error_message as $message)
		{
			echo "<li>{$message}</li>";
		}
		echo "</ul>
		</p>";
	}
	else
	{
		echo "<p>{$message}</p>";
	}

}
else
{
echo "
	<p>
		Use the form below to register a new account. Please be aware of that an activation e-mail will be sent when the account has been successfully registered. You need to take actions listed in the mail to be able to use your account.<br />
	</p>
	<form method=\"post\" action=\"./?p=register\" id=\"form_object\">
		<div>
			<ul>
				<li>
					<label for=\"email_1\">
						E-mail address
					</label>
					<input type=\"text\" name=\"email_1\" id=\"email_1\" />
				</li>
				<li>
					<label for=\"email_2\">
						Repeat e-mail address
					</label>
					<input type=\"text\" name=\"email_2\" id=\"email_2\" />
				</li>
				<li>
					<label for=\"pass_1\">
						Password
					</label>
					<input type=\"password\" name=\"pass_1\" id=\"pass_1\" />
				</li>
				<li>
					<label for=\"pass_2\">
						Repeat password
					</label>
					<input type=\"password\" name=\"pass_2\" id=\"pass_2\" />
				</li>
				<li class=\"secret\">
					<label for=\"foo\">
						Foo
					</label>
					<input type=\"text\" name=\"foo\" id=\"foo\" />
				</li>
			</ul>
			<ul>
				<li><input type=\"submit\" name=\"register\" value=\"Register\" /></li>
			</ul>
		</div>
	</form>";
}?>