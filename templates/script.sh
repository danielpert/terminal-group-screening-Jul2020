set -e
set -u

cd {{ project.config.project_dir }}

{% if group %}
{% for operation_group in operations|batch(group) %}
{% for operation in operation_group %}
{{ operation.cmd }} &
{% endfor %}
wait
{% endfor %}
{% else %}
{% for operation in operations %}
{{operation.cmd}}
{% endfor %}
{% endif %}
