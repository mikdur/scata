<?php
/**
 * Fördefinierade konstanter för meddelanden och returvärden
 */

/**
 * Sidspecifik information
 */
define("SITE_URL", "http://scata.mykopat.slu.se/");

/**
 * Sökvägar som filer ska sparas till
 */
define("DIR_DATASET_FAS", "/mykopat/scata/datasets");
define("DIR_DATASET_QUAL", "/mykopat/scata/datasets");
define("DIR_REFSEQ_FAS", "/mykopat/scata/referencesets");
define("DIR_TAGSET_TXT", "/mykopat/scata/tagsets");
define("DIR_JOBS", "/mykopat/scata/results");
define("DIR_OTHER", "/mykopat/scata/tmp");
define("DIR_TMP", "/mykopat/scata/tmp");

/**
 * Returvärden för registrering
 */
define("REGISTRATION_SUCCESS", 0);
define("REGISTRATION_NOT_MATCHING_PASSWORDS", 1);
define("REGISTRATION_NOT_VALID_PASSWORD", 2);
define("REGISTRATION_EMAIL_EXISTS", 3);
define("REGISTRATION_NOT_MATCHING_EMAILS", 4);
define("REGISTRATION_INVALID_ACTIVATION_KEY", 5);
define("REGISTRATION_EMAIL_DOESNT_EXIST",6);
define("REGISTRATION_MISSING_EMAIL", 7);
define("REGISTRATION_MISSING_UID", 8);

/**
 * Värden för användare
 */
define("USER_SUCCESS", 9);
define("USER_INVALID_CREDENTIALS", 10);
define("USER_NOT_ENABLED", 11);

/**
 * Värden för jobb
 */
define("JOB_SUCCESS", 12);
define("JOB_FAILURE", 13);

/**
 * Värden för filer
 */
define("FILE_DATASET_FAS", 14);
define("FILE_DATASET_QUAL", 15);
define("FILE_REFSET_FAS", 16);
define("FILE_TAGSET_TXT", 17);

/**
 * Värden för dataset
 */
define("DATASET_SUCCESS", 18);
define("DATASET_FAILURE", 19);



define("PASSWORD_SALT", "ai0sauli");
?>
