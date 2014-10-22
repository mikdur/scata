<?php
/**
 * Databashantering
 * @package TASP
 */

/**
 * Klass för databashantering.
 *
 * Använd {@link get_instance()} för att hämta ett databasobjekt. "new dbObject()" fungerar inte.
 *
 * Använd {@link query()} för att göra en fråga och {@link execute()} för att köra INSERT/UPDATE.<br>
 * {@link safe()} används för att skydda sig mot SQL-injections på införd data.<br>
 * Resultat hämtas med {@link fetch_data()}, {@link fetch_assoc()} eller {@link fetch_row()}.
 *
 * Klassen är definierad som en Singleton, dvs att endast ett databasobjekt kan skapas. Om objektet
 * redan är skapat när man kör {@link get_instance()} så returneras bara det befintliga objektet.
 *
 * Ändringar:
 *
 *	- 070808: Dokumentation /A.L.
 *	- 080624: Lägger till möjlighet att ange anslutningsuppgifter
 *			  i constructorn /A.L.
 *	- 090129: Tredje argumentet i safe() talar om ifall man ska köra addslashes() eller ej /A.L
 *	- 090620: Lånar {@link prepare()} från Wordpress 2.8 /A.L.
 *  - 090704: Konstanter för host, databas, användarnamn och lösenord hämtas från config.php.
 *			  dbObject omskriven som en Singleton. /A.L.
 *  - 100214: Tog bort konstaneter för host, databas, användarnamn och lösenord.
 * 
 * @example example_db.php Se exempel här
 * @uses config.php
 * @author Anders Lisspers (anders@lisspers.se)
 * @version 1.0
 * @package TASP
 */
class dbObject
{
	
/* Variabler
***************************************************************/
	private $host = "scata.mykopat.slu.se";
	private $database = "scata_dev";
	private $username = "scata_dev";
	private $password = "S4ta03v//";

	private $conn;
	private $debug;
	private $queries;
	
	/**
	 * Instansen av denna singleton
	 * @var string
	 */
	private static $the_instance;

	/**
     * Resultat från den senaste SQL-frågan
     * @var resource
     */
	public $results;
	
	/**
     * Antal rader som påverkades av senaste SQL-frågan
     * @var int
     */
	public $row_count;
	
	/**
     * Den senaste SQL-frågan
     * @var string
     */
	public $last_query;
	
	/**
     * ID för den rad som senast skapades med en INSERT-fråga
     * @var resource
     */
	public $last_id;


/* Konstruktor/destruktor
***************************************************************/
	/**
	 * Skapar anslutning till databasen
	 * @param string $host Adress till databasserver [optional]
	 * @param string $username Användarnamn [optional]
	 * @param string $pass Lösenord [optional]
	 * @param string $db Databasens namn [optional]
	 * @return void
	 */
	private function __construct($host = null, $username = null, $pass = null, $db = null)
	{
		$this->connect($host, $username, $pass, $db);
	}

	/**
	 * Stänger databasanslutningen.
	 * Om debug mode är på dumpas databasobjektet
	 * @return void
	 */
	public function __destruct()
	{
		if ($this->debug)
		{
			echo '<pre>';
			var_dump($this);
			echo '</pre>';
		}
	}

/* Metoder
***************************************************************/
	/**
	 * Skapar ett databasobjekt om det inte redan finns.
	 * Om ett objekt redan finns så returneras det
	 *
	 * @param string $host 
	 * @param string $username 
	 * @param string $pass 
	 * @param string $db 
	 * @return dbObject
	 */
	public static function get_instance($host = null, $username = null, $pass = null, $db = null)
	{
		if (!isset(self::$the_instance))
		{
			self::$the_instance = new dbObject($host, $username, $pass, $db);
		}
		
		return self::$the_instance;
	}

	/**
	 * Ansluter till angiven host och databas
	 * @return void
	 **/
	public function connect($host, $username, $password, $database)
	{
		if (isset($host) && isset($username) && isset($pass) && isset($db))
		{
			$this->host = $host;
			$this->username = $username;
			$this->password = $pass;
			$this->database = $db;
		}
		$this->conn = mysql_connect($this->host, $this->username, $this->password)
			or trigger_error('MySQL-error: ' . mysql_error(), E_USER_ERROR);
	
		mysql_select_db($this->database, $this->conn)
			or trigger_error('MySQL-error: ' . mysql_error(), E_USER_ERROR);
	
		mysql_query("SET time_zone = 'Europe/Stockholm'");
		
		
		$this->queries = 0;
	}
	
	/**
	 * Byter till en annan databas (på samma server)
	 * @param string $database Namn på den nya databasen
	 * @return void
	 */
	public function use_database($database)
	{
		mysql_select_db($database, $this->conn);
	}

	/**
	 * Utför en SQL-fråga (SELECT)
	 * @param string $sql Aktuell SQL-fråga
	 * @return resource
	 */
	public function query($sql)
	{
		if ($this->debug)
			echo "SQL-query: $sql<br />";

		$this->last_query = $sql;
		if (!is_resource($this->conn))
		{
			$this->set_debug(true);
			die();
		}
			
		$this->results = mysql_query($sql, $this->conn)
			or trigger_error('MySQL-error: ' . mysql_error() . "<br />Query: $sql", E_USER_ERROR);
		$this->row_count = mysql_num_rows($this->results);

		$this->queries++;
		
		return $this->results;
	}

	/**
	 * Utför en SQL-fråga (INSERT / UPDATE / DELETE)
	 * @param string $sql Aktuell SQL-fråga
	 * @return bool
	 */
	public function execute($sql)
	{
		if ($this->debug)
			echo "SQL-execution: $sql<br />";;
		
		$this->last_query = $sql;
		
		$this->results = mysql_query($sql, $this->conn)
			or trigger_error('MySQL-error: ' . mysql_error() . "<br />Query: $sql", E_USER_WARNING);
			
		$this->last_id = mysql_insert_id($this->conn);
		$this->queries++;
		$this->row_count = mysql_affected_rows();

		return $this->results;
	}
	
	/**
	 * Hämtar ut en rad ur ett frågeresultat och lägger i en numerisk array.
	 * @param resource $results  Resultat från en SQL-fråga
	 * @param int $row Om en specifik rad önskas anges den här
	 * @return array
	 */
	public function fetch_row($results = null, $row = -1)
	{
		if (!isset($results)) $results = $this->results;
		if ($row == -1) return mysql_fetch_row($results);
		for ($i = 0; $i <= $row; $i++)
		{
			$r = mysql_fetch_row($results);
		}
		return $r;
	}
	
	/**
	 * Hämtar ut en rad ur ett frågeresultat och lägger i en associativ array.
	 * @param resource $results  Resultat från en SQL-fråga
	 * @param int $row Om en specifik rad önskas anges den här
	 * @return array
	 */
	public function fetch_assoc($results = null, $row = -1)
	{
		if (!isset($results)) $results = $this->results;
		if ($row == -1) return mysql_fetch_assoc($results);
		for ($i = 0; $i <= $row; $i++)
		{
			$r = mysql_fetch_assoc($results);
		}
		return $r;
	}
	
	/**
	 * Hämtar ut data ur en SQL-fråga.
	 * Endast det första fältets värde returneras
	 * @param resource $results  Resultat från en SQL-fråga
	 * @return mixed
	 */
	public function fetch_data($results = null)
	{
		if (!isset($results)) $results = $this->results;
		if ($this->row_count > 0)
			return mysql_result($results, 0, 0);
		else
			return null;
	}
	

	/**
	 * Gör ett värde säkert att skicka i en SQL-fråga
	 * @param mixed $val Värdet som ska kontrolleras
	 * @param bool $add_slashes true om addslashes() ska köras oavsett
	 * @return mixed
	 */
	public function safe($val, $is_int = false, $add_slashes = false)
	{
		if ($is_int)
		{
			$id = intval($val);
			
			if ($id < 0)
				trigger_error("dbObject: Invalid input (id = $id < 0)", E_USER_ERROR);
				
			return $id;
		}
		else
		{		
		  return ($add_slashes || get_magic_quotes_gpc() != 1) ? addslashes($val) : $val;
		}
	}

	/**
	 * Function borrowed from Wordpress 2.8
     * Prepares a SQL query for safe execution.  Uses sprintf()-like syntax.
     *
     * This function only supports a small subset of the sprintf syntax; it only supports %d (decimal number), %s (string).
     * Does not support sign, padding, alignment, width or precision specifiers.
     * Does not support argument numbering/swapping.
     *
     * May be called like {@link http://php.net/sprintf sprintf()} or like {@link http://php.net/vsprintf vsprintf()}.
     *
     * Both %d and %s should be left unquoted in the query string.
     *
     * <code>
     * wpdb::prepare( "SELECT * FROM `table` WHERE `column` = %s AND `field` = %d", "foo", 1337 )
     * </code>
     *
     * @link http://php.net/sprintf Description of syntax.
     * @since 2.3.0
     *
     * @param string $query Query statement with sprintf()-like placeholders
     * @param array|mixed $args The array of variables to substitute into the query's placeholders if being called like {@link http://php.net/vsprintf vsprintf()}, or the first variable to substitute into the query's placeholders if being called like {@link http://php.net/sprintf sprintf()}.
     * @param mixed $args,... further variables to substitute into the query's placeholders if being called like {@link http://php.net/sprintf sprintf()}.
     * @return null|string Sanitized query string
     */
 	public function prepare($query = null)
	{
         if (is_null($query))
             return;
         $args = func_get_args();
         array_shift($args);

         // If args were passed as an array (as in vsprintf), move them up
         if (isset($args[0]) && is_array($args[0]))
             $args = $args[0];

         $query = str_replace("'%s'", '%s', $query); // in case someone mistakenly already singlequoted it
         $query = str_replace('"%s"', '%s', $query); // doublequote unquoting
         $query = str_replace('%s', "'%s'", $query); // quote the strings
         array_walk($args, array(&$this, 'safe'));
         return @vsprintf($query, $args);
     }


	/**
	 * Skriver ut innehållet i argumentet/argumenten.
	 * Innehållet skrivs ut med print_r, varefter exit() körs.
	 * Kan vara bra att använda i debugsyfte.
	 * dump() kan ta ett godtyckligt antal parametrar.
	 *
	 * @param mixed $dump,... Variabler som ska skrivas ut
	 * @return void
	 */
	private function dump()
	{
	    echo '<pre>';

	    $num_args = func_num_args();
	    for ($i = 0; $i < $num_args; ++$i)
	    {
	        print_r(func_get_arg($i));
	        echo "\n\n";
	    }

	    echo '</pre>';
	    exit;
	}
	

/* Setters
***************************************************************/
	/**
	 * Sätter / stänger av debug mode
	 * @param bool $debug 
	 * @return void
	 */
	public function set_debug($debug = true)
	{
		$this->debug = $debug;
	}
}
?>
