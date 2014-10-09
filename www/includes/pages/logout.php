<?php

$user = new User();
$user->logout();
redirect("./?p=home");
?>