<?php

/**
 * Klass för hantering av jobb
 */

class Job
{
	public $id;
	private $db;
	private $user;
	private $name;
	private $description;
	private $datasets;
	private $refseqs;
	private $paramset;
	private $filename;
	
	public function __construct($id = 0)
	{
		$this->db = dbObject::get_instance();
		$this->user = new User();
		$this->id = $id;

	}
	
	/**
	 * Anger id för job
	 * @param int id
	 * @return void
	 */
	public function set_id($id)
	{
		$this->id = $id;
	}
	
	/**
	 * Anger relaterade dataset
	 * @param array Datasets id
	 * @return void
	 */
	public function set_data_sets($datasets)
	{
		$this->datasets = $datasets;
	}
	
	/**
	 * Anger referenssekvenser
	 * @param array Referenssekvensers id
	 * @return void
	 */
	public function set_reference_sequences($refseqs)
	{
		$this->refseqs = $refseqs;
	}
	
	
	/**
	 * Sätter namn
	 * @param string Namn
	 * @return void
	 */
	public function set_name($name)
	{
		$this->name = $name;
	}
	
	/**
	 * Hämtar namn
	 * @return string Namn
	 */
	public function get_name()
	{
		if ($this->id == 0)
		{
			throw new Exception("Job: get_name(): id not set");
		}
		if (!isset($this->name))
		{
			$query = "SELECT name FROM Jobs WHERE idJobs = {$this->id}";
			$this->db->query($query);
			$this->name = $this->db->fetch_data();
		}
		return $this->name;
	}
	
	/**
	 * Sätter filnamn
	 * Anges inget filnamn blir namnet en slumpad alfanumeriskt sträng på 32 tecken
	 * Filnamnet är unikt
	 * @return void
	 */
	public function set_filename()
	{
		$success = false;
		do
		{
			$this->filename = rand_string();
			
			// Se om filnamn redan finns
			$query = "SELECT COUNT(idJobs) FROM Jobs WHERE filename = '{$this->filename}'";
			$this->db->query($query);
			if ($this->db->fetch_data() == 0)
			{
				$success == true;
			}
		}
		while(!$success);
	}
	
	
	/**
	 * Hämtar filnamn
	 * @return string Filnamn
	 */
	public function get_filename()
	{
		if ($this->id == 0)
		{
			throw new Exception("Job: get_filename(): id not set");
		}
		if (!isset($this->filename))
		{
			$query = "SELECT filename FROM Jobs WHERE idJobs = {$this->id}";
			$this->db->query($query);
			$this->filename = $this->db->fetch_data();
		}
		return $this->filename;
	}
	
	/**
	 * Sätter beskrivning
	 * @param string Beskrivning
	 * @return void
	 */
	public function set_description($description)
	{
		$this->description = $description;
	}
	
	/**
	 * Hämtar beskrivning
	 * @return string Beskrivning
	 */
	public function get_description()
	{
		if ($this->id == 0)
		{
			throw new Exception("Job: get_description(): id not set");
		}
		if (!isset($this->description))
		{
			$query = "SELECT description FROM Jobs WHERE idJobs = {$this->id}";
			$this->db->query($query);
			$this->description = $this->db->fetch_data();
		}
		return $this->description;
	}
		
	/**
	 * Sätter parameter set
	 * @param int ID för parameterset
	 * @return void
	 */
	public function set_parameter_set($paramset)
	{
		$this->paramset = $paramset;
	}
	
	/**
	 * Skapar ett job
	 * @return int JOB_SUCCESS om skapat, JOB_FAILURE annars
	 */
	public function create()
	{
		if (!isset($this->name))
		{
			throw new Exception("Job: create(): Name not set");
		}
		if (!isset($this->datasets) || !is_array($this->datasets))
		{
			throw new Exception("Job: create(): Dataset id(s) not set. Must be array.");
		}
		if (!isset($this->paramset))
		{
			throw new Exception("Job: create(): Parameter set id not set.");
		}
		if (!isset($this->description))
		{
			throw new Exception("Job: create(): Description not set.");
		}
		// Skapa jobb
		$query = "INSERT INTO
			Jobs(
				createTime,
				owner,
				status,
				jobSettingsSet,
				name,
				description)
			VALUES(
				NOW(),
				".$this->user->get_id().",
				1,
				{$this->paramset},
				'{$this->name}',
				'{$this->description}'
			)";
		$this->db->execute($query);
		$this->set_id($this->db->last_id);
		
		// Skapa relation jobb<->referens-set
		if (isset($this->refseqs) && is_array($this->refseqs))
		{
			for ($i = 0; $i < count($this->refseqs); $i++)
			{
				$query = "INSERT INTO
					JobRefs(
						refId,
						jobId)
					VALUES(
						{$this->refseqs[$i]},
						{$this->id}
						)";
				$this->db->execute($query);
			}
		}
		
		// Skapa relation jobb<->datast
		for ($i = 0; $i < count($this->datasets); $i++)
		{
			$query = "INSERT INTO
				JobDatasets(
					jobId,
					dataSetId)
				VALUES(
					{$this->id},
					{$this->datasets[$i]}
				)";
			$this->db->execute($query);
		}
		return JOB_SUCCESS;
	}
	
	/**
	 * Raderar ett job
	 * @return int True om raderat, false annars.
	 */
	public function delete()
	{
		if (!$this->is_owner() && !$this->user->is_admin())
		{
			return false;
		}

		// Radera relation jobb<->referens-sekvens
		$query = "DELETE FROM JobRefs WHERE jobId = {$this->id}";
		$this->db->execute($query);
		
		// Radera relation jobb<->dataset
		$query = "DELETE FROM JobDatasets WHERE jobId = {$this->id}";
		$this->db->execute($query);

		// Radera inlägg i job-tabellen
		$query = "DELETE FROM Jobs WHERE idJobs = {$this->id}";
		$this->db->execute($query);
		return true;
		
	}
	
	/**
	 * Kontrollerar att användare är ägare av jobbet
	 * @return boolean true användare är ägare, false annars
	 */
	public function is_owner()
	{
		if (!isset($this->id) || $this->id == 0)
		{
			throw new Exception("Job: is_owner(): Id not set.");
		}
		$this->db->query("SELECT owner = {$this->user->get_id()} FROM Jobs WHERE idJobs = {$this->id}");
		return ($this->db->fetch_data() == 1);
	}
	
	/**
	 * Ser om användare har rättighet till jobbet
	 * @return boolean True om användare har det, false annars.
	 */
	public function user_has_rights()
	{
		if (!isset($this->id) || $this->id == 0)
		{
			throw new Exception("Job: is_owner(): Id not set.");
		}
		return ($this->is_owner() || $this->user->is_admin());
	}
	
	/**
	 * Kontrollerar om jobbet är klart
	 * @return boolean true om klart, false annars
	 */
	public function is_ready()
	{
		if (!isset($this->id) || $this->id == 0)
		{
			throw new Exception("Job: is_owner(): Id not set.");
		}
		$this->db->query("SELECT status = 3 FROM Jobs WHERE idJobs = {$this->id}");
		return $this->db->fetch_data();
	}
	
	/**
	 * Hämtar status för jobb
	 * @return string Status uttryckt i text
	 */
	public function get_status()
	{
		if (!isset($this->id) || $this->id == 0)
		{
			throw new Exception("Job: get_status(): Id not set.");
		}
		$query = "SELECT
			JobStatus.Status
		FROM
			JobStatus
		INNER JOIN
			Jobs
			ON
				Jobs.status = JobStatus.idJobStatus
		WHERE
			Jobs.idJobs = {$this->id}";
		$this->db->query($query);
		return $this->db->fetch_data();
	}
	
	/**
	 * Hämtar id för associerade dataset
	 * @return array ID för datasets
	 */
	public function get_datasets()
	{
		if (!isset($this->id) || $this->id == 0)
		{
			throw new Exception("Job: get_datasets(): Id not set.");
		}
		$query = "SELECT dataSetId FROM JobDatasets WHERE jobId = {$this->id}";
		$this->db->query($query);
		
		$datasets = Array();
		while($dataset = $this->db->fetch_row())
		{
			$datasets[] = $dataset[0];
		}
		return $datasets;
	}
	
	/**
	 * Hämtar id för associerade referessekvenser
	 * @param array ID för referenssekvens
	 */
	public function get_refseqs()
	{
		if (!isset($this->id) || $this->id == 0)
		{
			throw new Exception("Job: get_refseqs(): Id not set.");
		}
		$query = "SELECT refId FROM JobRefs WHERE jobId = {$this->id}";
		$this->db->query($query);
		
		$refseqs = Array();
		while($refseq = $this->db->fetch_row())
		{
			$refseqs[] = $refseq[0];
		}
		return $refseqs;
	}
	
	/**
	 * Hämtar id parameterset
	 * @return int ID för settet
	 */
	public function get_parameterset()
	{
		if (!isset($this->id) || $this->id == 0)
		{
			throw new Exception("Job: get_status(): Id not set.");
		}
		$query = "SELECT jobSettingsSet FROM Jobs WHERE idJobs = {$this->id}";
		$this->db->query($query);
		return $this->db->fetch_data();
	}
	

}

?>
