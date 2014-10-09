<?php
/**
 * Klass för att generera formulär
 * @package TASP
 */

/**
 * Skapar och hanterar HTML-formulär, med validering och säkerhet.
 * Formulärsfält läggs till med {@link add_field()} (obs, läs beskrivning),
 * javascript för felhantering genereras med {@link generate_javascript()} och
 * html-koden för formuläret genereras med {@link generate_form()}.
 *
 * Ändringar:
 *	- 081222: Objekten får unika id (möjligt att ha flera FormObject på en sida) /A.L.
 *	- 081229: Döper om generate_javascript() och generate_form() till print_* /A.L
 *	- 081229: Lägger till set_action(), set_submit_label() och use_groups(). /A.L
 *	- 081230: Bara de två första argumenten är obligatoriska i {@link add_field()} /A.L.
 *  - 081230: Man kan nu använda must_be=X och print=false i {@link add_field} /A.L.
 *	- 090107: Lade till stöd för URL-validering /A.L.
 *	- 090509: Nytt felmeddelandesystem och stöd för selectionboxar /A.L.
 *	- 090701: * När nåt fel uppstått, och formuläret skrivs ut igen, markeras de felaktiga fälten med css-klassen "form_error".
 *	          * {@link validates()} kör {@link handle_input_data()} innan fälten valideras.
 *	          * Ny metod för antispam: {@link enable_antispam()}.
 *	          /A.L.
 *	- 090715: Stöd för multiselect i selectboxar, och stöd för method="get". /A.L.
 *	- 090716: Möjlighet att ange onclick="" och onchange="" som options till {@link add_field()} /A.L.
 *
 * @example example_form.php Se exempel här
 * @package TASP
 * @uses dbObject
 * @author Anders Lisspers <anders@lisspers.se>
 * @version 1.0
 */
class FormObject
{

/* Variabler
***************************************************************/
	private $db;
	private $db_table;
	private $form_id;
	private $form_fileupload;
	private $antispam_enabled;
	
	private $id_field;
	private $id_value;
	
	private $debug;
	
	private $use_groups = false;
	private $groups = array();
	public $fields = array();
	private $select_options = array();
	private $errors = array();
	private $error_fields = array();
	
	private $form_method = 'post';
	private $form_action = null;
	private $submit_label = null;
	
	// Namn på krypteringsfunktion för lösenord
	private $password_encryption = 'sha1';
	
	// Varje formulärsobjekt måste ha olika id
	private $obj_id;
	// Räknar upp varje gång en instans skapas
	public static $next_obj_id = 0;
	
	private $text_no = 0;

	
/* Konstruktor/destruktor
***************************************************************/
	/**
	 * Skapar en ny instans av FormObject.
	 * De två första argumenten är obligatoriska.
	 * @param string $table Databastabell som används för användarinformation
	 * @param bool $debug true om innehållet i $this ska dumpas vid destruct
	 * @return void
	 */
	public function __construct($table = null, $debug = false)
	{
		$this->db = dbObject::get_instance();
		
		if (isset($table))
			$this->db_table = $table;
		else
			trigger_error('FormObject: No database table given.', E_USER_ERROR);
		
		$this->obj_id = self::$next_obj_id++;
		
		$this->debug = $debug;
	}
	
	/**
	 * Funktion som körs när instansen tas bort.
	 * Om $this->debug är satt till true dumpas innehållet i instansen
	 * @return void
	 */
	public function __destruct()
	{
		if ($this->debug)
		{
			echo get_class($this) . ':<br /><pre>';
			var_dump($this);
			echo '</pre>';
		}
	}


/* Metoder
***************************************************************/
	/**
	 * Kontrollerar om formuläret har postats
	 * @return true om formuläret har postats, false annars
	 */
	public function got_post()
	{
		return isset($_POST['form_submit' . $this->obj_id]);
	}

	/**
	 * Fyller flera fältvärden samtidigt.
	 *
	 * Två alternativ för anrop till denna funktion:
	 * 1) Argumentet är en associativ array där nyckeln är fältets namn i databasen och
	 *    värdet är dess värde ($arr['fältets_namn'] = 'värde').
	 *    Arrayen loopas igenom och värden läggs in i respektive fält.
	 * 2) Inget argument skickas med. Då måste {@link id_field} och {@link id_value}} vara satt för den rad
	 *    som är tänkt att fylla formuläret.
	 * 
	 * @param array|string $arr Se de olika alternativet ovan
	 * @param string $id_value Se ovan [Optional]
	 * @return void
	 */
	public function fill_values($arr = null)
	{
		if (is_array($arr))
		{
			# Alternativ 1: Fyll formuläret med angiven data
			foreach ($arr as $key => $value)
			{
				foreach ($this->fields as $id => $field)
				{
					if ($field['field_name'] == $key)
						$this->set_value($id, $value);
				}
			}
		}
		else
		{
			# Alternativ 2: Hämta data från databasen och fyll formuläret med denna
			if (!isset($this->id_field) || !isset($this->id_value))
				trigger_error('FormObject (fill_values()): ID field and ID value must be set (run set_id_field/value()-methods).', E_USER_ERROR);
			
			$db = $this->db;
			
			$db->query("SELECT {$this->get_db_fields()} FROM {$this->db_table} WHERE {$this->id_field} = '{$this->id_value}'");
			if ($db->row_count != 1)
			{
				return;
			}
			else
			{
				$this->fill_values($db->fetch_assoc());
			}
			
		}
	}

	/**
	 * Gör formuläret redo till att ladda upp filer
	 * @return void
	 */
	public function enable_fileupload()
	{
		$this->form_fileupload = true;
	}

	/**
	 * Lägger till ett textfält som döljs via css. Fältet är tänkt att inte fyllas i.
	 * Spamrobotar fyller alltid i alla fält som finns i ett formulär, så om detta (dolda)
	 * fält är ifyllt måste således formuläret fyllts i av en robot.
	 *
	 * Idé tagen från http://swedishfika.com/2009/05/29/are-you-tired-of-captcha/
	 *
	 * @return void
	 **/
	public function enable_antispam()
	{
		$this->antispam_enabled = true;
	}

	/**
	 * Tar emot införda värden
	 * Sparar infört värde till respektive fält med {@link set_value()}
	 * @return void
	 */
	public function handle_input_data()
	{
		$db = $this->db;
		if ($this->antispam_enabled)
		{
			if (strlen($_POST['form_foolspam']) != 0)
			{
				trigger_error('FormObject: The system believes you\'re a spamrobot. If you\'re not, please contact the administrator.', E_USER_ERROR);
			}
		}
		
		foreach ($this->fields as $id => $field)
		{
			if ($field['type'] == 'just_text')
			{
				continue;
			}
			
			# Erik! Detta �r det som strular. F�ltet blir aldrig satt till en array n�r man vill det.
			if (isset($field['options']['mandatory'])){
			  if(!isset($_POST[$id])){
			    $this->fields[$id]['value'] = 'MANDATORY_NOT_SET';
			    continue;
			  }elseif(is_array($_POST[$id]) && 
				  count($_POST[$id]) == 0){
			    $this->fields[$id]['value'] = 'MANDATORY_NOT_SET';
			    continue;
			  }elseif(!is_array($_POST[$id]) && 
				  strlen($_POST[$id]) == 0){
			    $this->fields[$id]['value'] = 'MANDATORY_NOT_SET';
			    continue;
			  }
			}

			if (!isset($_POST[$id]))
			{
				if ($this->fields[$id]['type'] == 'checkbox')
				{
					$this->fields[$id]['value'] = '';
				}
				else
				{
					continue;
				}
			}
			elseif (is_object($this->db))
			{
			  if(is_array($_POST[$id])){
			    foreach ($_POST[$id] as $key => $value) { 
			      $new_array[$key] = addslashes($_POST[$id][$key]); 
			    }
			    $this->fields[$id]['value'] = $new_array;
			  }else{
			    $this->fields[$id]['value'] = $db->safe($_POST[$id]);
			  }
			}
			else
			{
				$this->fields[$id]['value'] = $_POST[$id];

			}

		}
	}
	
	/**
	 * Lägger in en ny fältgrupp
	 * @param string id Unikt namn på gruppen
	 * @param string text Texten som ska stå som rubrik för gruppen
	 * @return void
	 */
	public function add_group($id, $text = '')
	{
		$this->groups[$id] = array('text' => $text, 'members' => array());
	}
	
	/**
	 * Lägger in en textrad i formuläret
	 * @param string text Texten som ska läggas till
	 * @return void
	 */
	public function add_text($text)
	{
		$this->add_field('text' . $this->text_no++, '', $text, 'just_text');
	}
	
	/**
	 * Lägger till ett alternativ till en selectbox
	 * @param string $id Namnet på select-fältet
	 * @param string $value Värde på alternativet (dvs <option value="$value">)
	 * @param string $text Text som ska stå för alternativet (dvs <option>$text</option>)
	 * @param boolean $selected true om alternativet ska vara förvalt [Optional]
	 * @return void
	 */
	public function add_option($id, $value, $text, $selected = false) {
		$this->select_options[$id][$value]['text'] = $text;
		$this->select_options[$id][$value]['selected'] = $selected;
	}
	
	/**
	 * Lägger till ett formulärsfält
	 *
	 * ************************************************************************************<br>
	 * BESKRIVNING AV INSTÄLLNINGAR FÖR FORMULÄRSFÄLTEN<br>
	 *
	 * _______$type______________________________________________________________________<br>
	 * Det finns åtta typer av fält man kan lägga till. Sju av dessa är vanliga formulärstyper:
	 * text, textarea, checkbox, radio, select, password och file.
	 *
	 * Observera att för checkbox krävs "checkbox_value" i options, för radio krävs "radio_alternatives"
	 * och "radio_values" och för select krävs att man kör {@link add_option()}:
	 *
	 * Den sjunde används om man vill ha annan html-kod än dessa.
	 * Skriv då in html-koden istället för type.
	 *
	 * _______$options___________________________________________________________________<br>
	 * När man lägger till ett fält (med ld) finns en parameter som heter $options.
	 * Denna ska vara en semikolonseparerad lista med olika flaggor/inställningar.
	 *
	 * Följande inställningar finns tillgängliga:<br>
	 *	- minlength=X			-	Tvingar att indatan är minst X tecken lång.
	 *	- unique=true			-	Angivet värde får inte finnas i databasen sedan tidigare.
	 *	- mandatory=true		-	Fältet är obligatoriskt.
	 * 								Skriver ut '<span class="form_mandatory">*</span>' efter fälttexten
	 *	- print=false			-	Om fältet inte ska skrivas ut
	 *	- label_after=X			-	Text som ska stå bakom formulärsfältet
	 *	- confirm=true			-	Två textrutor ska fyllas i och matcha (confirm_label krävs).
	 *	- confirm_label=X		-	Sätter texten vid bekräfta-textrutan till X.
	 *	- checkbox_value=X		-	Sätter värdet som ska skickas med, om checkboxen är vald, till X
	 *	- radio_alternatives=X	-	| (pipe)-separerad lista med radiobutton-alternativ
	 *	- radio_values=X		-	| (pipe)-separerad lista med värden på radiobutton-alternativen
	 *	- multiselect=true		-	Används för selectboxar där man ska kunna välja flera alternativ
	 *	- class=X				-	CSS-klass. Lägger till class="X" i html-tagen.
	 *	- must_be=X				-	Värdet på fältet måste vara lika med X
	 *	- onclick=X				-	Sätter onclick="X" i fältet
	 *	- onchange=X			-	Sätter onchange="X" i fältet
	 *	- validate_as=X			-	Vilken typ av validering som ska användas.<br>
	 *
	 *	För tillgängliga valideringstyper, se {@link is_valid()}.
	 *	
	 * Exempel:
	 * 	add_field('username','u_username','Användarnamn: ','text','minlength=6;confirm=true;confirm_label=Bekräfta:', 'grupp-id');
	 *
	 * **********************************************************************************
	 *
	 *
	 * @param string id Unikt namn på fältet
	 * @param string field_name Namn på fältet i databasen
	 * @param string label Text som beskriver fältet (t.ex. "Användarnamn: ")
	 * @param string type Fältets typ (se ovan).
	 * @param string options Se options-beskrivning ovan
	 * @return void
	 */
	public function add_field($id, $field_name, $label = null, $type = null, $options = '', $group = null)
	{
		if (isset($this->fields[$id]))
			trigger_error("FormObject: Given ID is already occupied ($id)", E_USER_ERROR);
			
		if (strlen($options))
		{
			$opt_array = explode(';', $options);
			$new_opt_array = array();
			foreach ($opt_array as $opt)
			{
				$opt_name = substr($opt, 0, strpos($opt, '='));
				$opt_value = substr($opt, strpos($opt, '=') + 1);
				$new_opt_array[$opt_name] = $this->fix_str($opt_value); 
			}
		}
		else
		{
			$new_opt_array = null;
		}
		$opts = $new_opt_array;
		
		
		if (isset($opts['confirm']))
		{
			if (!isset($opts['confirm_label']))
				trigger_error('FormObject: Confirm requires confirm_label in $options (' . $id . ')', E_USER_ERROR);
			$confirm_label = $opts['confirm_label'];
			unset($opts['confirm']);
			unset($opts['confirm_label']);
						
			$add_confirm = true;
		}
		

		$this->fields[$id] = array('field_name' => $field_name,
								   'label' => $label,
								   'type' => $type,
								   'options' => $opts
							 );
		if (isset($group))
			$this->groups[$group]['members'][] = $id;
					
		if (isset($add_confirm) && $add_confirm)
		{
			if (isset($opts['mandatory']))
			{
				$this->add_field($id . '_confirm', '', $this->fix_str($confirm_label), $type, 'mandatory=true;equal_to=' . $id, $group);
			}
			else
			{
				$this->add_field($id . '_confirm', '', $this->fix_str($confirm_label), $type, 'equal_to=' . $id, $group);
			}
		}
	}
	
	/**
	 * Tar bort ett formulärsfält
	 * @param string id Id på fältet som ska tas bort
	 * @return void
	 */
	public function remove_field($id)
	{
		unset($this->fields[$id]);
	}
	

	/**
	 * Skriver ut html-kod för formuläret
	 * @param string action Dit formulärdatan ska postas [optional]
	 * @param string submit_label Text på submitknappen [optional]
	 * @param bool use_groups Om true sorteras fälten efter dess grupper [optional]
	 * @return void
	 */
	public function print_form($action = null, $submit_label = null, $use_groups = false)
	{
		if (!$use_groups && $this->use_groups)
			$use_groups = true;
			
		echo $this->generate_form($action, $submit_label, $use_groups);
	}
	
	/**
	 * Genererar html-kod för formuläret.
	 * @param string action Dit formulärdatan ska postas [optional]
	 * @param string submit_label Text på submitknappen [optional]
	 * @param bool use_groups Om true sorteras fälten efter dess grupper [optional]
	 * @return string HTML-kod för formulär
	 */
	public function generate_form($action = null, $submit_label = null, $use_groups = false)
	{
		if ($action == null)
		{
			if ($this->form_action == null)
				trigger_error('FormObject (generate_form()): Form action not set.', E_USER_ERROR);
			else
				$action = $this->form_action;
		}	
		if ($submit_label == null)
		{
			if ($this->submit_label == null)
				trigger_error('FormObject (generate_form()): Submit label not set.', E_USER_ERROR);
			else
				$submit_label = $this->submit_label;
		}
		
		$form_id = (isset($this->form_id)) ? $this->form_id : 'form_object';
		$ret_val = "
		<form id=\"$form_id\" method=\"" . $this->form_method . "\" action=\"$action\"" . ($this->form_fileupload ? ' enctype="multipart/form-data"' : '') . " onsubmit=\"return validate_form$this->obj_id();\">";
		
		if ($use_groups)
		{
			foreach ($this->groups as $id => $group)
			{
				  
				$ret_val .= "
					<fieldset class=\"form_fieldset\" id=\"form_grp_$id\">";
					
				if (isset($group['text']) && strlen($group['text']) > 0)
					$ret_val .= "
						{$group['text']}";

				$ret_val .= '<ul>';
			
				foreach ($group['members'] as $field_id)
				{
					$ret_val .= $this->generate_field($field_id, $this->fields[$field_id]);
				}
				$ret_val .= '
			</ul>
			</fieldset>';
			}
		}
		else
		{
			$ret_val .= '
			<ul>';
			foreach ($this->fields as $id => $field)
			{
				$ret_val .= $this->generate_field($id, $field);
			}
			$ret_val .= '
			</ul>';
		}
		
		$ret_val .= '
			<ul>
				<li>
					<input type="submit" id="form_submit' . $this->obj_id . '" name="form_submit' . $this->obj_id . '" value="' . $submit_label . '" />
					' . ($this->antispam_enabled ? '<label class="form_hide" for="foolspam">Don\'t write anything here!</label><input type="text" class="form_hide" id="form_foolspam" name="form_foolspam" value="" />' : '') . '
				</li>
			</ul>
		</form>';
		
		return $ret_val;
	}
	
	
	/**
	 * Genererar htmlkod för aktuellt formulärsfält
	 * @param string $id Namn på aktuellt fält
	 * @param array $field Aktuellt fält
	 * @return void
	 */
	private function generate_field($id, $field)
	{
		unset($options);
		
		$error_field = in_array($id, $this->error_fields);
		
		if (!isset($ret_val))
			$ret_val = '';
			
		if (isset($field['options']))
			$options = $field['options'];
		
		if ($field['type'] != 'hidden' && (!isset($options['print']) || $options['print'] != 'false'))
			$ret_val .= '
			<li>';

		if (strlen($field['label']) > 0 && $field['type'] != 'checkbox' && $field['type'] != 'just_text')
		{
			$ret_val .= '
				<label' . ($error_field ? ' class="form_error"' : '') . ' id="form_' . $id . '_label" for="form_' . $id . '">
					' . $field['label'] .
					(isset($options['mandatory']) ? ' <span class="form_mandatory">*</span>' : '') . '
				</label>';
		}
		
		$parameters = ' id="form_' . $id . '" name="' . $id . '"' .
			(isset($options['class']) ? ' class="' . $options['class'] . '"' : '') .
			(isset($options['onchange']) ? ' onchange="' . $options['onchange'] . '"' : '') . 
			(isset($options['onclick']) ? ' onclick="' . $options['onclick'] . '"' : '') .
			(isset($options['default']) ? ' value="' . $options['default'] . '"' : '');

		switch ($field['type'])
		{
			case 'text':
			case 'password':
			case 'file':
				$ret_val .= '
				<input' . $parameters . ' type="' . $field['type'] . '" ';
				
			 	if (isset($options['minlength']) || isset($options['equal_to']))
				{
					$ret_val .= 'onkeyup="';

					if (isset($options['minlength']))
						$ret_val .= "check_length('$id', " . $options['minlength'] . ");";

					if (isset($options['equal_to']))
						$ret_val .= "check_equal('$id', '" . $options['equal_to'] . "');";

					$ret_val .= '"';
				}

				if (isset($field['value']) && $field['type'] != 'password')
					$ret_val .= ' value="' . $field['value'] . '"';

				if (isset($options['class']))
					$ret_val .= ' class="' . $options['class'] . '"';

				$ret_val .= ' />';

				if (isset($options['label_after']))
					$ret_val .= $options['label_after'];
				break;

			case 'checkbox':
				if (!isset($options['checkbox_value']))
					trigger_error('FormObject: Checkbox requires checkbox_value=X in $options. (' . $id . ')', E_USER_ERROR);

				$ret_val .= '
				<label class="form_normal' . ($error_field ? ' form_error' : '') . '" id="form_' . $id . '_label" for="form_' . $id . '">
					<input' . $parameters . ' type="checkbox" value="' . $options['checkbox_value'] . '"' . (isset($field['value']) && $field['value'] == $options['checkbox_value'] ? ' checked="checked"' : '') . ' />
					' . $field['label'] . '</label>' . 
					(isset($options['mandatory']) && $options['mandatory'] == 'true' ? ' <span class="form_mandatory">*</span>' : '');
				break;
			
			case 'radio':
				if (!isset($options['radio_alternatives']) || !isset($options['radio_values']))
					trigger_error('FormObject: Radio requires radio_alternatives=X and radio_values=Y in $options (' . $id . ')', E_USER_ERROR);
				
				$alternatives = explode('|', $options['radio_alternatives']);
				$values = explode('|', $options['radio_values']);
				
				if (count($alternatives) != count($values))
					trigger_error('FormObject: Number of radio_alternatives must be equal to number of radio_values (' . $id . ')', E_USER_ERROR);
					
				for ($i = 0; $i < count($alternatives); $i++)
				{
					$ret_val .= "
					<label class=\"form_normal" . ($error_field ? ' form_error' : '') . "\" for=\"form_{$id}_{$i}\" id=\"form_{$id}_{$i}_label\">
						<input id=\"form_{$id}_{$i}\" name=\"$id\" type=\"radio\" value=\"$values[$i]\"" . (isset($field['value']) && $field['value'] == $values[$i] ? ' checked="checked"' : '') . " />
						$alternatives[$i]
					</label>";
				}
				
				break;
				
			case 'textarea':
				$ret_val .= '
				<textarea' . $parameters . ' rows="" cols="" onkeyup="';
				if (isset($options['minlength']))
					$ret_val .= "check_length('$id', " . $options['minlength'] . ");";
				$ret_val .= '">';
				if (isset($field['value']))
					$ret_val .= $field['value'];					
				$ret_val .= '</textarea>';
				if (isset($options['label_after']))
					$ret_val .= $options['label_after'];
				break;
			
			case 'select':
				if (!isset($this->select_options[$id]) || !is_array($this->select_options[$id]))
				{
					trigger_error("FormObject: No options set for selectbox ($id)", E_USER_ERROR);
				}
					
				$ret_val .= "
				<select" .
					(isset($options['multiselect']) ? str_replace("name=\"{$id}\"", "name=\"{$id}[]\"", $parameters . ' multiple="multiple"') : $parameters) .
				">";
				foreach ($this->select_options[$id] as $value => $option)
				{
					$ret_val .= "
					<option value=\"{$value}\"" . ($option['selected'] || (isset($field['value']) && $field['value'] == $value) ? ' selected="selected"' : '') . ">{$option['text']}</option>";
				}
				$ret_val .= "
				</select>";
				
				if (isset($options['label_after']))
					$ret_val .= $options['label_after'];
				break;
				
			case 'just_text':
				$ret_val .= $field['label'];
				break;

			case 'hidden':
				$ret_val .= '<input ' . $parameters . ' type="hidden"' . (isset($field['value']) ? ' value="' . $field['value'] . '"' : '') . ' />';
				break;

			default:
				$ret_val .= '
				' . $field['type'];

				if (isset($options['label_after']))
					$ret_val .= $options['label_after'];					
				break;
		}

		if ($field['type'] != 'hidden')
			$ret_val .= '
			</li>';
		
		return $ret_val;
	}
	
	/**
	 * Skriver ut kod för javascriptvalidering
	 * @return string
	 **/
	public function print_javascript()
	{	
		echo $this->generate_javascript();
	}
	
	/**
	 * Genererar kod för javascriptvalidering
	 * @return string
	 **/
	public function generate_javascript()
	{	
		$ret_val = "
		<script type=\"text/javascript\">
		//<![CDATA[	
		function validate_form{$this->obj_id}()
		{
			if (!document.getElementById) return true;
			var form_error = false;
			var filter = [];
			filter['email'] = /^[_a-zA-Z0-9-]+(\.[_a-zA-Z0-9-]+)*@[a-zA-Z0-9-]+(\.[a-zA-Z0-9-]+)*(\.[a-z]{2,6})$/;
			filter['url'] = /^(https?):\/\/((?:[-a-z0-9]+\.)*)([-a-z0-9]+)\.(com|edu|gov|int|mil|net|org|biz|info|name|museum|coop|aero|[a-z]{2})(\/[-a-z0-9_:@&?=+,.!\/~*%$]*)?$/;
			filter['phone'] = /^[\+(?:\d{1,3}|0]+[\)\(\d\s-]*?$/;
			filter['number'] = /^([0-9]+)$/;
			filter['date'] = /^(19|20)?\d{2}[-|\s]?(0[0-9]|1[0-2])[-|\s]?([0-2][0-9]|3[0|1])$/;
		";

		foreach ($this->fields as $id => $field)
		{
			unset($options);
			if (isset($field['options']))
				$options = $field['options'];
			
			if (isset($options['minlength']) || isset($options['mandatory']) || isset($options['validate_as']) || isset($options['equal_to']))
			{
				$statements = '';

				if (isset($options['minlength']))
					$statements .= 'elem.value.length < ' . $options['minlength'];
					
				if (isset($options['mandatory']))
				{
					switch($field['type'])
					{
						case 'text':
						case 'password':
						case 'textarea':
						case 'file':
							$statements .= (strlen($statements) ? ' || ' : '') . 'elem.value.length == 0';
							break;
						case 'checkbox':
							$statements .= (strlen($statements) ? ' || ' : '') . '!elem.checked';
							break;
					}
				}
				
				if (isset($options['equal_to']))
				{
					$statements .= (strlen($statements) > 0 ? ' || ' : '') . "!check_equal('{$id}', '{$options['equal_to']}')";
				}
				
				if (isset($options['validate_as']))
				{
					if (!isset($options['mandatory']))
					{
						$statements .= (strlen($statements) > 0 ? ' || ' : '') . 'elem.value.length > 0 && !filter[\'' . $options['validate_as'] . '\'].test(elem.value)';
					}
					else
					{
						$statements .= (strlen($statements) > 0 ? ' || ' : '') . '!filter[\'' . $options['validate_as'] . '\'].test(elem.value)';
					}
				}
				
				if (strlen($statements) > 0)
				{
					$ret_val .= "
				// Validerar $id
					var elem = document.getElementById('form_$id');
					if ($statements)
					{
						document.getElementById('form_" . $id . "_label').className = 'form_error';
						form_error = true;
					}
					else
					{
						document.getElementById('form_" . $id . "_label').className = '';
					}
					";
				}
			}
		}	
			
		$ret_val .= '
			if (form_error)
			{
				alert("An error occurred. Please correct the marked fields.");
				return false;
			}
			else return true;
		}
		
		function check_length(id, length)
		{
			id = "form_" + id;
			if (document.getElementById(id).value.length < length)
				document.getElementById(id + "_label").className = "form_error";
			else
				document.getElementById(id + "_label").className = "";
		}
		
		function check_equal(id1, id2)
		{
			id1 = "form_" + id1;
			id2 = "form_" + id2;
			var elem1 = document.getElementById(id1);
			var elem2 = document.getElementById(id2);
			var elem1label = document.getElementById(id1 + "_label");
			var elem2label = document.getElementById(id2 + "_label");
			
			if (elem1.value != elem2.value)
			{
				elem1label.className = "form_error";
				return false;
			}
			else
			{
				elem1label.className = "";
				return true;
			}
		}
		//]]>	
		</script>';
		return $ret_val;
	}
	
	/**
	 * Alias till funktionen {@link validate_form()}, dock körs {@link handle_input_data()} innan.
	 * @return bool
	 */
	public function validates()
	{
		$this->handle_input_data();
		
		return $this->validate_form();
	}
	
	/**
	 * Validerar införd data.
	 * Data valideras enligt de options som valts för aktuellt fält.
	 * Om unique=true görs en sql-fråga mot databasen för att se om
	 * värdet finns innan eller ej.
	 *
	 * Det som kontrolleras är:
	 * 	'validation' 	- T.ex. att ett telefonnummer är på rätt form.
	 * 	'must_be'		- Att värdet stämmer med det värde som det ska vara.
	 *	'minlength'		- Att värdet har nått minimumlängden
	 *	'mandatory'		- Har det obligatoriska fältet något värde?
	 *	'unique'		- Är värdet unikt i databasen?
	 *	'equal_to'		- Är värdet lika som det andra fältet?
	 *
	 * Om någon kontroll inte validerar, t.ex
	 * om värdet på fältet "epost" inte är en giltig epostadress, läggs
	 * 'epost_validation' in i arrayen $this->errors.
	 * Om ett fält "msg" inte är tillräckligt långt läggs
	 * 'msg_minlength' in i arrayen.
	 *
	 * Uppstår något fel returnerar validate_form() false, och då
	 * kan felen hämtas genom {@link get_errors()} eller skrivas ut
	 * med {@link print_errors()}.
	 *
	 * @return bool true om all data validerar, false annars
	 */
	public function validate_form()
	{
		$db = $this->db;
		$valid = true; 
		foreach ($this->fields as $id => $field)
		{	
			$this_valid = true;
			unset($options);
			
			if (isset($field['options']['validate_as']))
			{
				if (isset($field['options']['mandatory']))
				{
					$this_valid = $this->is_valid($field['value'], $field['options']['validate_as']);
				}
				else
				{
					$this_valid = strlen($field['value']) == 0 || $this->is_valid($field['value'], $field['options']['validate_as']);
				}
			}	
			if (!$this_valid)
			{
				$valid = false;
				$this->error_fields[] = $id;
				$this->errors[] = "{$id}_validation";
				if ($field['value'] == 'MANDATORY_NOT_SET')
					$this->set_value($id, '');
			}
			elseif (isset($field['options']))
			{
				$options = $field['options'];
				
				foreach ($options as $option => $value)
				{
					switch ($option) 
					{
						case 'must_be':
							if ($field['value'] != $value)
							{
								$this->error_fields[] = $id;
								$this->errors[] = "{$id}_must_be";
								$valid = false;
								continue;
							}
							break;

						case 'minlength':
							if (strlen($field['value']) < $value)
							{
								$this->error_fields[] = $id;
								$this->errors[] = "{$id}_minlength";
								$valid = false;
								continue;
							}
							break;
							
						case 'mandatory':
						  if ($field['value'] == 'MANDATORY_NOT_SET')
							{
								$this->error_fields[] = $id;
								$this->errors[] = "{$id}_mandatory";
								$this->set_value($id, '');
								$valid = false;
								continue;
							}
							break;
							
						case 'unique':
							$sql = $db->safe("SELECT " . $field['field_name'] . " FROM $this->db_table WHERE " . $field['field_name'] . " = \"" . $field['value'] . "\"");
							$db->query($sql);
							if ($db->row_count)
							{
								$this->error_fields[] = $id;
								$this->errors[] = "{$id}_unique";
								$valid = false;
								
							}
							if (!$valid) continue;
							break;
							
						case 'equal_to':
							if ($field['value'] != $this->fields[$value]['value'])
							{
								$this->error_fields[] = $id;
								$this->errors[] = "{$id}_confirm";
								$valid = false;
								continue;
							}
							
							break;
					}
				}
			}
		}
		return $valid;
	}
	
	
	/**
	 * Kontrollerar om värdet från ett formulärsfält är giltigt.
	 * Tillgängliga typer av kontrollering:
	 *	- 'phone' 	- Telefonnummer
	 *	- 'pers_nr' - Personnummer
	 *	- 'email' 	- E-postadress
	 *	- 'url'		- URL
	 *	- 'number' 	- Siffra
	 *  - 'date'	- Datum
	 * @param mixed $val Värdet som ska kontrolleras
	 * @param string $type Vilken typ av kontroll som ska ske
	 * @return bool
	 */
	private function is_valid($val, $type)
	{
		switch (strtolower($type))
		{
			case 'text':
				return true;
				break;
				
			case 'phone':
				return preg_match('/^[\+(?:\d{1,3}|0]+[\)\(\d\s-]*?$/', $val);
				break;

			case 'pers_nr':
				$match = preg_match('/^(?:19|20)?\d{6}-?\d{4}$/', $val);
				if (!$match) return false;
				$val = preg_replace('/^(?:19|20)?(\d{6})-?(\d{4})/', '\\1\\2', $val);
				$n = 2;
				$sum = 0;

				for ($i=0; $i < 9; $i++)
				{
					$tmp = $val[$i] * $n;
					($tmp > 9) ? $sum += 1 + ($tmp % 10) : $sum += $tmp;
					($n == 2) ? $n = 1 : $n = 2;
				}

				return !(($sum + $val[9]) % 10);
				break;

			case 'email':
				return preg_match('/^[_a-zA-Z0-9-]+(\.[_a-zA-Z0-9-]+)*@[a-zA-Z0-9-]+(\.[a-zA-Z0-9-]+)*(\.[a-z]{2,6})$/', $val);
				break;
			
			case 'url':
				return preg_match('/^(https?):\/\/((?:[-a-z0-9]+\.)*)([-a-z0-9]+)\.(com|edu|gov|int|mil|net|org|biz|info|name|museum|coop|aero|[a-z]{2})(\/[-a-z0-9_:@&?=+,.!\/~*%$]*)?$/i', $val);
				break;

			case 'number':
				$match = preg_match('/^([0-9]+)$/', $val);
				if (isset($min_length) && $min_length)
					return ($match && strlen($val) >= $min_length);
				return $match;
				break;
			
			case 'date':                                                                       
				return preg_match('/^(19|20)?\d{2}[-|\s]?(0[0-9]|1[0-2])[-|\s]?([0-2][0-9]|3[0|1])$/', $val);
				break;
		}
		return false;
	}
	
	/**
	 * Skriver ut felmeddelanden som uppstod vid validering.
	 * Felen skrivs ut som en lista (ul) med css-klassen "form_error"
	 *
	 * Egna meddelanden skickas med som en array till funktionen
	 * där alla nycklar motsvarar något fel som "skapas" i
	 * {@link validate_form()}.
	 *
	 * Ex: Fältet med id "epost" hade inte ett värde som är
	 * en giltig e-postadress. Därför sparades felet
	 * "epost_validation" vid validate_form(). I argumentet
	 * till denna funktion bör då ett värde finnas som har
	 * nyckeln "epost_validation". Dvs till exempel:
	 * $my_error_msgs["epost_validation"] = "Du angav en felaktig e-postadress";
	 *
	 * Om ett fel hittas som inte har något angivet felmeddelande i argumentet
	 * skrivs meddelandet "Uncaught error (epost_validation)" ut.
	 *
	 * @param $err_msgs Arrayen med egna felmeddelanden
	 * @return void
	 */
	public function print_errors($err_msgs = null)
	{
		if (count($this->errors) == 0)
			return;
			
		echo '<ul class="form_error">';
		foreach ($this->errors as $error)
		{
			echo '<li>';

			if (isset($err_msgs[$error]))
				echo $err_msgs[$error];
			else
				echo "Uncaught error ({$error})";

			echo '</li>';
		}
		echo '</ul>';
	}
	
	
	/**
	 * Sparar införd data i en ny rad i databasen
	 * @return void
	 */
	public function insert()
	{
		$db_fields = '';
		$db_values = '';
		
		foreach ($this->fields as $field)
		{
			if ($field['type'] != 'just_text')
			{	
				$db_fields .= $field['field_name'] . ', ';
				$value = $field['value'];
				if ($field['type'] == 'password' && strlen($this->password_encryption = ''))
				{
					$encrypt = $this->password_encryption;
					if (!function_exists($encrypt))
					{
						trigger_error('FormObject (insert()): Password encryption function not found ("' . $encrypt . '").', E_USER_ERROR);
					}
					
					$value = $encrypt($value);
				}
				$db_values .= "'" . $field['value'] . "', ";
			}
		}
		$db_fields = substr($db_fields, 0, strlen($db_fields) - 2);
		$db_values = substr($db_values, 0, strlen($db_values) - 2);
		$sql = "INSERT INTO $this->db_table ($db_fields) VALUES ($db_values)";
		$this->db->execute($sql);
	}
	
	/**
	 * Uppdaterar en befintlig rad i databasen med införd data
	 * @param string $id_field Fält i databasen som identifierar rätt rad
	 * @param mixed $id_value Värdet på id-fältet för den rad som ska uppdateras
	 * @return void
	 */
	public function update($id_field = null, $id_value = null)
	{
		if (!isset($id_field))
		{
			$id_field = $this->id_field;
		}
		
		if (!isset($id_value))
		{
			$id_value = $this->id_value;
		}
		
		$db_values = '';
		foreach ($this->fields as $field)
		{
			if ($field['type'] != 'just_text')
			{
				$value = $field['value'];
				if ($field['type'] == 'password' && strlen($this->password_encryption) > 0)
				{
					$encrypt = $this->password_encryption;
					if (!function_exists($encrypt))
					{
						trigger_error('FormObject (insert()): Password encryption function not found ("' . $encrypt . '").', E_USER_ERROR);
					}
					
					$value = $encrypt($value);
				}
				$db_values .= $field['field_name'] . ' = \'' . $value . '\', ';
				
				
			}
		}
		$db_values = substr($db_values, 0, strlen($db_values) - 2);
		$sql = "UPDATE $this->db_table SET $db_values WHERE $id_field = '" . $this->db->safe($id_value) . "'";
		$this->db->execute($sql);
	}
	
	/**
	 * Fixar svenska tecken till html-entities
	 * @param string $str Sträng att fixa
	 * @return string
	 */
	private function fix_str($str)
	{
		$str = str_replace('å', 'å', $str);
		$str = str_replace('ä', 'ä', $str);
		$str = str_replace('ö', 'ö', $str);
		
		return $str;		
	}
	

/* Getters
***************************************************************/
	/**
	 * Hämtar värdet från ett fält
	 * @return string
	 */
	public function get_value($id)
	{
		return $this->fields[$id]['value'];
	}

	/**
	 * Hämtar ut felen som uppstod vid validering
	 * @return array
	 */
	public function get_errors()
	{
		return $this->errors;
	}
	

	/**
	 * Hämtar ut id på senast införda rad i databasen
	 * @return int
	 */
	public function get_insert_id()
	{
		return $this->db->last_id;
	}

	/**
	 * Hämtar ut databasfält för givet fält
	 * @param string $field Fält
	 * @return string
	 */
	public function get_db_field($field)
	{
		return ($this->fields[$field]['field_name']);
	}

	/**
	 * Hämtar ut alla databasfält i kommaseparerad lista (användbar till SQL-frågor)
	 * @return string
	 */
	public function get_db_fields()
	{
		foreach ($this->fields as $id => $field)
		{
			if ($field['type'] != 'just_text')
				$fields[] = $this->get_db_field($id);
		}
		
		return implode(', ', $fields);
	}

/* Setters
***************************************************************/
	/**
	 * Sätter css-klassen på ett fält
	 * @param string $id Fält-id
	 * @param string $class Klassnamn
	 */
	public function set_class($id, $class)
	{
		$this->fields[$id]['options']['class'] = $class;
	}
	
	/**
	 * Sätter värdet för ett fält
	 * @param string $id Fält-id
	 * @param string $str Värde
	 * @return void
	 */
	public function set_value($id, $str)
	{
		$this->fields[$id]['value'] = $str;
	}
	
	/**
	 * Sätter metod på formuläret
	 * @param string $method "post" eller "get"
	 * @return void
	 */
	public function set_method($method)
	{
		$this->form_method = $method;
	}
	
	/**
	 * Sätter action på formuläret
	 * @param string $action URL
	 * @return void
	 */
	public function set_action($action)
	{
		$this->form_action = $action;
	}

	/**
	 * Sätter texten på submitknappen
	 * @return void
	 */
	public function set_submit_label($submit)
	{
		$this->submit_label = $submit;
	}

	/**
	 * Sätter om grupper ska användas eller ej
	 * @param bool $groups
	 * @return void
	 */
	public function use_groups($groups = true)
	{
		$this->use_groups = $groups;
	}

	/**
	 * Sätter css-id på formuläret
	 * @return void
	 */
	public function set_form_id($id)
	{
		$this->form_id = $id;
	}
	
	/**
	 * Sätter fält i databasen som identifierar aktuell rad
	 * @return void
	 */
	public function set_id_field($field)
	{
		$this->id_field = $this->db->safe($field);
	}
	
	/**
	 * Sätter värdet på det fält i databasen som identifierar aktuell rad
	 * @return void
	 */
	public function set_id_value($value)
	{
		$this->id_value = $this->db->safe($value);
	}
	
	/**
	 * Sätter vilken krypteringsfunktion som ska användas till lösenord
	 * @param string $func Namn på krypteringsfunktion
	 * @return void
	 */
	public function set_password_encryption($func)
	{
		$this->password_encryption = $func;
	}
	
	/**
	 * Sätter $debug-variabeln till true eller false
	 * @param bool $debug
	 * @return void
	 */
	public function set_debug($debug = true)
	{
		$this->debug = $debug;
	}
	
	
}
?>