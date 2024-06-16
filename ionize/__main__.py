import click
from .Solution import Solution
from .Database import Database


@click.group()
def cli():
    pass


@cli.command()
def database():
    '''
    Print the entire ion database.
    '''
    click.echo(Database().serialize())


@cli.command()
@click.argument('name')
def ion(name):
    '''
    Print the physical properties of an ion.
    '''
    db = Database()
    if name in db:
        result = db[name]
        click.echo(result.serialize(nested = False, compact = True))
    elif (search_results := db.search(name)):
        click.echo('Did you mean one of these?')
        for result in search_results:
            click.echo(f"\t{result}")
    else:
        click.echo(f"'{name}' not found")
    


@cli.command()
@click.option('--component', '-c', type=(str, float), multiple=True)
@click.option('--titrate', '-t', type = (str, float), default = (None, None), help = 'titrate the solution with the specified titrant to the specified value of a property, e.g. \'hydrochloric acid\' 8.5')
@click.option('--titration_property', '-p', default = 'pH', help = 'physical property to titrate, default pH')
@click.option('--temperature', default = '25', help = "The temperature at which to simulate the solution, in Celcius.")
def solution(component, titrate, titration_property, temperature):
    '''
    Compute the physical properties of a solution.
    
    IONS: comma-separated list of ion names, e.g. 'tris,hydrochloric acid'
    CONCENTRATIONS: comma-separated list of concentrations (molar) in same order as ions, e.g. '0.01,0.006'
    '''
    
    ions = [ion for ion, concentration in component]
    concentrations = [concentration for ion, concentration in component]
    initial_solution = Solution(ions, concentrations)
    initial_solution.temperature(temperature)
    
    if titrate != (None, None): click.echo('starting solution:')
    click.echo(initial_solution)
    click.echo(initial_solution.serialize(nested = False, compact = True))
    
    if titrate != (None, None):
        titrant, target = titrate
        titrated_solution = initial_solution.titrate(titrant, target, titration_property)
        click.echo('\ntitrated solution:')
        click.echo(titrated_solution)
        click.echo(titrated_solution.serialize(nested = False, compact = True))
        click.echo('\ntitration: add %f M %s' % ((titrated_solution - initial_solution).concentration(titrant), titrant))


if __name__ == '__main__':
    cli()

