<?php
/**
 * Klass för hantering av dataset
 */

class Dataset
{
	public $id;
	private $db;
	private $user;
	private $date;
	private $owner;
	private $name;
	private $primer5;
	private $primer5match;
	private $primer3;
	private $primer3match;
	private $tagset5;
	private $tagset3;
	private $min_len;
	private $mean_qual;
	private $min_qual;
	private $max_len;
	private $file_type;
	private $raw_filtering;
	
	public function __construct($id = 0)
	{
		$this->db = dbObject::get_instance();
		$this->user = new User();
		$this->name;
		$this->id = $id;
	}
	
	/**
	 * Sätt namn datasettet
	 * @param name Namn
	 * @return void
	 */
	public function set_name($name)
	{
		$this->name = $name;
	}
	
	/**
	 * Hämtar namn på settet
	 * @return string Namn
	 */
	public function get_name()
	{
		if ($this->id == 0)
		{
			throw new Exception("Dataset: get_name(): id not set");
		}
		if (!isset($this->name))
		{
			$query = "SELECT Name FROM Datasets WHERE idDatasets = {$this->id}";
			$this->db->query($query);
			$this->name = $this->db->fetch_data();
		}
		return $this->name;
	}

	/**
	 * Sätt primer5 f�r dataset
	 * @param primer5 Primer5
	 * @return void
	 */
	public function set_primer5($primer5)
	{
		$this->primer5 = $primer5;
	}
	
	/**
	 * Hämtar primer5 på settet
	 * @return string primer5
	 */
	public function get_primer5()
	{
		if ($this->id == 0)
		{
			throw new Exception("Dataset: get_primer5(): id not set");
		}
		if (!isset($this->primer5))
		{
			$query = "SELECT Primer5 FROM Datasets WHERE idDatasets = {$this->id}";
			$this->db->query($query);
			$this->primer5 = $this->db->fetch_data();
		}
		return $this->primer5;
	}

	/**
	 * Sätt primer5score f�r dataset
	 * @param primer5score Primer5score
	 * @return void
	 */
	public function set_primer5score($primer5score)
	{
		$this->primer5score = $primer5score;
	}
	
	/**
	 * Hämtar primer5score på settet
	 * @return string primer5score
	 */
	public function get_primer5score()
	{
		if ($this->id == 0)
		{
			throw new Exception("Dataset: get_primer5score(): id not set");
		}
		if (!isset($this->primer5score))
		{
			$query = "SELECT Primer5score FROM Datasets WHERE idDatasets = {$this->id}";
			$this->db->query($query);
			$this->primer5score = $this->db->fetch_data();
		}
		return $this->primer5score;
	}

	/**
	 * Sätt primer5 f�r dataset
	 * @param primer3 Primer3
	 * @return void
	 */
	public function set_primer3($primer3)
	{
		$this->primer3 = $primer3;
	}
	
	/**
	 * Hämtar primer3 på settet
	 * @return string primer3
	 */
	public function get_primer3()
	{
		if ($this->id == 0)
		{
			throw new Exception("Dataset: get_primer3(): id not set");
		}
		if (!isset($this->primer3))
		{
			$query = "SELECT Primer3 FROM Datasets WHERE idDatasets = {$this->id}";
			$this->db->query($query);
			$this->primer3 = $this->db->fetch_data();
		}
		return $this->primer3;
	}

	/**
	 * Sätt primer3score f�r dataset
	 * @param primer3score Primer3score
	 * @return void
	 */
	public function set_primer3score($primer3score)
	{
		$this->primer3score = $primer3score;
	}
	
	/**
	 * Hämtar primer3score på settet
	 * @return string primer3score
	 */
	public function get_primer3score()
	{
		if ($this->id == 0)
		{
			throw new Exception("Dataset: get_primer3score(): id not set");
		}
		if (!isset($this->primer3score))
		{
			$query = "SELECT Primer3score FROM Datasets WHERE idDatasets = {$this->id}";
			$this->db->query($query);
			$this->primer3score = $this->db->fetch_data();
		}
		return $this->primer3score;
	}

	/**
	 * Sätt tagset5 f�r dataset
	 * @param tagset5 Tagset5
	 * @return void
	 */
	public function set_tagset5($tagset5)
	{
		$this->tagset5 = $tagset5;
	}
	
	/**
	 * Hämtar tagset5 på settet
	 * @return string tagset5
	 */
	public function get_tagset5()
	{
		if ($this->id == 0)
		{
			throw new Exception("Dataset: get_tagset5(): id not set");
		}
		if (!isset($this->tagset5))
		{
			$query = "SELECT Tagset5 FROM Datasets WHERE idDatasets = {$this->id}";
			$this->db->query($query);
			$this->tagset5 = $this->db->fetch_data();
		}
		return $this->tagset5;
	}

	/**
	 * Sätt tagset3 f�r dataset
	 * @param tagset3 Tagset3
	 * @return void
	 */
	public function set_tagset3($tagset3)
	{
		$this->tagset3 = $tagset3;
	}
	
	/**
	 * Hämtar tagset3 på settet
	 * @return string tagset3
	 */
	public function get_tagset3()
	{
		if ($this->id == 0)
		{
			throw new Exception("Dataset: get_tagset3(): id not set");
		}
		if (!isset($this->tagset3))
		{
			$query = "SELECT Tagset3 FROM Datasets WHERE idDatasets = {$this->id}";
			$this->db->query($query);
			$this->tagset3 = $this->db->fetch_data();
		}
		return $this->tagset3;
	}

	/**
	 * Sätt min_len f�r dataset
	 * @param min_len Min_Len
	 * @return void
	 */
	public function set_min_len($min_len)
	{
		$this->min_len = $min_len;
	}
	
	/**
	 * Hämtar min_len på settet
	 * @return string min_len
	 */
	public function get_min_len()
	{
		if ($this->id == 0)
		{
			throw new Exception("Dataset: get_min_Len(): id not set");
		}
		if (!isset($this->min_len))
		{
			$query = "SELECT min_len FROM Datasets WHERE idDatasets = {$this->id}";
			$this->db->query($query);
			$this->min_len = $this->db->fetch_data();
		}
		return $this->min_len;
	}

	/**
	 * Sätt max_len f�r dataset
	 * @param max_len Max_Len
	 * @return void
	 */
	public function set_max_len($max_len)
	{
		$this->max_len = $max_len;
	}
	
	/**
	 * Hämtar max_len på settet
	 * @return string max_len
	 */
	public function get_max_len()
	{
		if ($this->id == 0)
		{
			throw new Exception("Dataset: get_max_len(): id not set");
		}
		if (!isset($this->max_len))
		{
			$query = "SELECT max_len FROM Datasets WHERE idDatasets = {$this->id}";
			$this->db->query($query);
			$this->name = $this->db->fetch_data();
		}
		return $this->max_len;
	}

	/**
	 * Sätt mean_qual f�r dataset
	 * @param mean_qual Mean_Qual
	 * @return void
	 */
	public function set_mean_qual($mean_qual)
	{
		$this->mean_qual = $mean_qual;
	}
	
	/**
	 * Hämtar mean_qual på settet
	 * @return string mean_qual
	 */
	public function get_mean_qual()
	{
		if ($this->id == 0)
		{
			throw new Exception("Dataset: get_mean_qual(): id not set");
		}
		if (!isset($this->mean_qual))
		{
			$query = "SELECT mean_qual FROM Datasets WHERE idDatasets = {$this->id}";
			$this->db->query($query);
			$this->name = $this->db->fetch_data();
		}
		return $this->mean_qual;
	}

	/**
	 * Sätt min_qual f�r dataset
	 * @param min_qual Min_Qual
	 * @return void
	 */
	public function set_min_qual($min_qual)
	{
		$this->min_qual = $min_qual;
	}
	
	/**
	 * Hämtar min_qual på settet
	 * @return string min_qual
	 */
	public function get_min_qual()
	{
		if ($this->id == 0)
		{
			throw new Exception("Dataset: get_min_qual(): id not set");
		}
		if (!isset($this->min_qual))
		{
			$query = "SELECT min_qual FROM Datasets WHERE idDatasets = {$this->id}";
			$this->db->query($query);
			$this->name = $this->db->fetch_data();
		}
		return $this->min_qual;
	}
	
	/**
	 * Sets file type
	 */
	public function set_file_type($file_type)
        {
                $this->file_type = $file_type;
        }

	/**
	 * Gets file type
	 */
        public function get_file_type()
        {
                if ($this->id == 0)
                {
                        throw new Exception("Dataset: get_file_type(): id not set");
                }
                if (!isset($this->file_type))
                {
                        $query = "SELECT file_type FROM Datasets WHERE idDatasets = {$this->id}";
                        $this->db->query($query);
                        $this->file_type = $this->db->fetch_data();
                }
                return $this->file_type;
        }



	/**
	 * Sätt raw_filtering f�r dataset
	 * @param raw_filtering Raw_Filtering
	 * @return void
	 */
	public function set_raw_filtering($raw_filtering)
	{
		$this->raw_filtering = $raw_filtering;
	}
	
	/**
	 * Hämtar raw_filtering på settet
	 * @return string raw_filtering
	 */
	public function get_raw_filtering()
	{
		if ($this->id == 0)
		{
			throw new Exception("Dataset: get_raw_filtering(): id not set");
		}
		if (!isset($this->raw_filtering))
		{
			$query = "SELECT raw_filtering FROM Datasets WHERE idDatasets = {$this->id}";
			$this->db->query($query);
			$this->raw_filtering = $this->db->fetch_data();
		}
		return $this->raw_filtering;
	}


	/**
	 * Hämtar description på settet
	 * @return string min_qual
	 */
	public function get_description()
	{
		if ($this->id == 0)
		{
			throw new Exception("Dataset: get_description(): id not set");
		}
		if (!isset($this->description))
		{
			$query = "SELECT Description FROM Datasets WHERE idDatasets = {$this->id}";
			$this->db->query($query);
			$this->description = $this->db->fetch_data();
		}
		return $this->description;
	}

	
	/**
	 * Hämtar datum (YYYY-mm-dd) för skapandet av datasettet
	 * @return string Datum
	 */
	public function get_created_date()
	{
		if ($this->id == 0)
		{
			throw new Exception("Dataset: get_created_date(): id not set");
		}
		if (!isset($this->date))
		{
			$query = "SELECT DATE(createdDate) FROM Datasets WHERE idDatasets = {$this->id}";
			$this->db->query($query);
			$this->date = $this->db->fetch_data();
		}
		return $this->date;
	}
	
	/**
	 * Skapa dataset
	 * Sätter $this->id till datasettets id
	 * @return int DATASET_SUCCESS om skrivningen skedde. DATASET_FAILURE om användare ej inloggad.
	 */
	public function create()
	{
		if ($this->user->get_id() == 0)
		{
			throw new Exception("Dataset: create: user id not set");
		}
		if (!isset($this->name) || $this->name == "")
		{
			throw new Exception("Dataset: create(): Name not set");
		}
		$query = "INSERT INTO Datasets(
			createdDate,
			Name,
			Owner,
			ready,
			Primer5,
			Primer3,
			Primer5score,
			Primer3score,
			Tagset5,
			Tagset3,
			min_len,
			max_len,
			mean_qual,
			min_qual,
                        raw_filtering,
			file_type)
		VALUES(
			NOW(),
			'{$this->name}',
			".$this->user->get_id().",
			0,
			'{$this->primer5}',
			'{$this->primer3}',
			'{$this->primer5score}',
			'{$this->primer3score}',
			'{$this->tagset5}',
			'{$this->tagset3}',
			'{$this->min_len}',
			'{$this->max_len}',
			'{$this->mean_qual}',
			'{$this->min_qual}',
                        '{$this->raw_filtering}',
			'{$this->file_type}')";
		$this->db->execute($query);
		
		$this->id = $this->db->last_id;
		
		return DATASET_SUCCESS;
	}

	/**
	 * Ser om dataset är redo
	 * @return boolean True om redo, false annars.
	 */
	public function is_ready()
	{
		if ($this->id == 0)
		{
			throw new Exception("Dataset: is_ready(): id not set");
		}
			$query = "SELECT ready > 0 FROM Datasets WHERE idDatasets = {$this->id}";
		$this->db->query($query);
		return $this->db->fetch_data();
	}

	/**
	 * Ser om användare är ägare
	 * @return boolean True om ägare, false annars.
	 */
	public function is_owner()
	{
		if ($this->id == 0)
		{
			throw new Exception("Dataset: is_owner(): id not set");
		}
		if (!isset($this->owner))
		{
			$query = "SELECT owner FROM Datasets WHERE idDatasets = {$this->id}";
			$this->db->query($query);
			$this->owner = $this->db->fetch_data();
		}
		return $this->owner == $this->user->get_id();
	}

	/**
	 * Ta bort dataset
	 * Anledningar till att dataset inte kan tas bort är:
	 * 1. Användaren är inte admin och inte ägare av settet
	 * 2. Settet är kopplat till ett jobb
	 * @param int id
	 * @return boolean True om raderad, false annars.
	 */
	public function delete($id = 0)
	{
		if ($id == 0)
		{
			$id = $this->id;
		}
		if ($id == 0)
		{
			throw new Exception("Dataset: delete(): Id not set");
		}
		// Se om användare har rättigheter att radera
		$query = "SELECT ".($this->user->is_admin() ? 1 : 0)." = 1 OR owner = ".$this->user->get_id()." FROM Datasets WHERE idDatasets = {$this->id}";
		$this->db->query($query);
		if ($this->db->fetch_data() == 0)
		{
			echo "No rights";
			return false;
		}
		
		// Se om settet är kopplat till ett jobb
		$query = "SELECT COUNT(idJobDatset) > 0 FROM JobDatasets WHERE dataSetId = {$this->id}";
		$this->db->query($query);
		if ($this->db->fetch_data() == 1)
		{
			echo "Relation exists";
			return false;
		}

		$query = "DELETE FROM Datasets WHERE idDatasets = {$this->id}";
		
		// Radera filen
		return $this->db->execute($query) && @unlink(DIR_DATASET_QUAL ."/". $this->id . ".seqs.pick") && @unlink(DIR_DATASET_FAS ."/". $this->id . ".stat.pick");
	}
}
?>
