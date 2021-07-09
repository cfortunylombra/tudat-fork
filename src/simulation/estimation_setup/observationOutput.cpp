/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/simulation/estimation_setup/observationOutput.h"

namespace tudat
{

namespace simulation_setup
{


void checkObservationDependentVariableEnvironment(
        const SystemOfBodies& bodies,
        const std::shared_ptr< StationAngleObservationDependentVariableSettings > variableSettings )
{
    if( isObservationDependentVariableGroundStationProperty( variableSettings ) )
    {
//        if( variableSettings->relevantLinkEnds_.size( ) != 1 )
//        {
//            throw std::runtime_error( "Error in observation dependent variables when creating function for " +
//                                      getObservationDependentVariableName( variableSettings ) +
//                                      ", only one input link end can be processed" );
//        }

        std::string bodyName = variableSettings->relevantLinkEnd_.first;
        std::string stationName = variableSettings->relevantLinkEnd_.second;

        if( bodies.count( bodyName ) == 0 )
        {
            throw std::runtime_error( "Error in observation dependent variables when creating function for " +
                                      getObservationDependentVariableName( variableSettings ) +
                                      ", could not find body " + bodyName );
        }

        if( bodies.at( bodyName )->getGroundStationMap( ).count( stationName ) == 0 )
        {
            throw std::runtime_error( "Error in observation dependent variables when creating function for " +
                                      getObservationDependentVariableName( variableSettings ) +
                                      ", could not find station " + stationName + " on body " + bodyName );
        }
    }
}

double getLinkEndRange(
        const std::vector< Eigen::Matrix< double, 6, 1 > >& linkEndStates,
        const std::pair< int, int > linkEndIndices )
{
    return ( linkEndStates.at( linkEndIndices.first ).segment( 0, 3 ) -
                linkEndStates.at( linkEndIndices.second ).segment( 0, 3 ) ).norm( );
}

DoubleObservationDependentVariableFunction getBodyAvoidanceFunction(
        const SystemOfBodies& bodies,
        const std::shared_ptr< BodyAvoidanceObservationDependentVariableSettings > variableSettings,
        const observation_models::ObservableType observableType,
        const observation_models::LinkEnds linkEnds )
{
    // Check if observation is integrated
    bool isObservableIntegrated = observation_models::isObservableOfIntegratedType( observableType );
    if( isObservableIntegrated == 1 )
    {
        throw std::runtime_error( "Error in body-avoidance observation dependent variables for " +
                                  getObservationDependentVariableName( variableSettings ) +
                                  " Interface for integrated observables not yet implemented" );
    }

    int numberOfLinks = observation_models::getNumberOfLinksInObservable(
                observableType, linkEnds.size( ) );

    std::pair< int, int > linkEndIndices;
    DoubleObservationDependentVariableFunction outputFunction;
    if( numberOfLinks == 1 )
    {
        linkEndIndices = std::make_pair( 0, 1 );
    }


    outputFunction  = [=]( const std::vector< double >& linkEndTimes,
            const std::vector< Eigen::Matrix< double, 6, 1 > >& linkEndStates,
            const Eigen::VectorXd& observationValue )
    {
        double averageTime = ( linkEndTimes.at( linkEndIndices.first ) +
                linkEndTimes.at( linkEndIndices.second ) ) / 2.0;
        Eigen::Vector3d bodyToAvoid =
                bodies.at( variableSettings->bodyAvoidance_ )->getStateInBaseFrameFromEphemeris< double, double >( averageTime ).segment< 3 >( 0 );
        return observation_models::computeCosineBodyAvoidanceAngle( linkEndStates, linkEndIndices, bodyToAvoid );
    };
    return outputFunction;
}

DoubleObservationDependentVariableFunction getTargetRangeFunction(
        const std::shared_ptr< InterlinkObservationDependentVariableSettings > variableSettings,
        const observation_models::ObservableType observableType,
        const observation_models::LinkEnds linkEnds )
{
    // Check if observation is integrated
    bool isObservableIntegrated = observation_models::isObservableOfIntegratedType( observableType );
    if( isObservableIntegrated == 1 )
    {
        throw std::runtime_error( "Error in target-range observation dependent variables for " +
                                  getObservationDependentVariableName( variableSettings ) +
                                  " Interface for integrated observables not yet implemented" );
    }

    int numberOfLinks = observation_models::getNumberOfLinksInObservable(
                observableType, linkEnds.size( ) );

    std::pair< int, int > linkEndIndices;
    DoubleObservationDependentVariableFunction outputFunction;
    if( numberOfLinks == 1 )
    {
        linkEndIndices = std::make_pair( 0, 1 );
    }

    outputFunction  = [=]( const std::vector< double >&,
            const std::vector< Eigen::Matrix< double, 6, 1 > >& linkEndStates,
            const Eigen::VectorXd& observationValue )
    { return getLinkEndRange( linkEndStates, linkEndIndices ); };
    return outputFunction;
}

DoubleObservationDependentVariableFunction getStationObservationAngleFunction(
        const SystemOfBodies& bodies,
        const std::shared_ptr< StationAngleObservationDependentVariableSettings > variableSettings,
        const observation_models::ObservableType observableType,
        const observation_models::LinkEnds linkEnds )
{
    checkObservationDependentVariableEnvironment( bodies, variableSettings );

    // Retrieve link-end ID for station
    std::string bodyName = variableSettings->relevantLinkEnd_.first;
    std::string stationName = variableSettings->relevantLinkEnd_.second;

    // Check if observation is integrated
    bool isObservableIntegrated = observation_models::isObservableOfIntegratedType( observableType );
    if( isObservableIntegrated == 1 )
    {
        throw std::runtime_error( "Error in station-angle observation dependent variables for " +
                                  getObservationDependentVariableName( variableSettings ) +
                                  " Interface for integrated observables not yet implemented" );
    }

    // Get link-end indices consistent with settings
    std::vector< std::pair< int, int > > stateTimeIndex =
            observation_models::getLinkStateAndTimeIndicesForLinkEnd(
                linkEnds, observableType, variableSettings->relevantLinkEnd_ );
    if( stateTimeIndex.size( ) == 0 )
    {
        throw std::runtime_error( "Error in station-angle observation dependent variables for " +
                                  getObservationDependentVariableName( variableSettings ) +
                                  " Could not find link end: (" + variableSettings->relevantLinkEnd_.first + ", " +
                                  variableSettings->relevantLinkEnd_.second + ")" );
    }

    std::pair< int, int > linkEndIndicesToUse;
    if( stateTimeIndex.size( ) == 1 )
    {
        linkEndIndicesToUse = stateTimeIndex.at( 0 );
    }
    else if( stateTimeIndex.size( ) > 1 )
    {
        throw std::runtime_error( "Error in station-angle observation dependent variables for " +
                                  getObservationDependentVariableName( variableSettings ) +
                                  " Interface not yet implemented for multiple viable link indices" );
    }


    // Retrieve pointing angles calculator for station, and setup output function
    std::shared_ptr< ground_stations::PointingAnglesCalculator > pointingAnglesCalculator =
            bodies.at( bodyName )->getGroundStationMap( ).at( stationName )->getPointingAnglesCalculator( );
    DoubleObservationDependentVariableFunction outputFunction;
    if( variableSettings->variableType_ == station_elevation_angle )
    {
        outputFunction = [=]( const std::vector< double >& linkEndTimes,
                const std::vector< Eigen::Matrix< double, 6, 1 > >& linkEndStates,
                const Eigen::VectorXd& observationValue )
        { return ground_stations::calculateGroundStationElevationAngle(
                        pointingAnglesCalculator, linkEndStates, linkEndTimes, linkEndIndicesToUse ); };
    }
    else if( variableSettings->variableType_ == station_azimuth_angle )
    {
        outputFunction = [=]( const std::vector< double >& linkEndTimes,
                const std::vector< Eigen::Matrix< double, 6, 1 > >& linkEndStates,
                const Eigen::VectorXd& observationValue )
        { return ground_stations::calculateGroundStationAzimuthAngle(
                        pointingAnglesCalculator, linkEndStates, linkEndTimes, linkEndIndicesToUse ); };
    }
    return outputFunction;
}


DoubleObservationDependentVariableFunction getObservationDoubleDependentVariableFunction(
        const SystemOfBodies& bodies,
        const std::shared_ptr< ObservationDependentVariableSettings > variableSettings,
        const observation_models::ObservableType observableType,
        const observation_models::LinkEnds linkEnds )
{

    DoubleObservationDependentVariableFunction outputFunction;
    switch( variableSettings->variableType_ )
    {
    case station_elevation_angle:
    {
        std::shared_ptr< StationAngleObservationDependentVariableSettings > angleSettings =
                std::dynamic_pointer_cast< StationAngleObservationDependentVariableSettings >(
                    variableSettings );
        if( angleSettings == nullptr )
        {
            throw std::runtime_error( "Error in observation dependent variables, incorrect type found for station_elevation_angle" );
        }
        outputFunction = getStationObservationAngleFunction(
                    bodies, angleSettings, observableType, linkEnds );
        break;
    }
    case station_azimuth_angle:
    {
        std::shared_ptr< StationAngleObservationDependentVariableSettings > angleSettings =
                std::dynamic_pointer_cast< StationAngleObservationDependentVariableSettings >(
                    variableSettings );
        if( angleSettings == nullptr )
        {
            throw std::runtime_error( "Error in observation dependent variables, incorrect type found for station_azimuth_angle" );
        }

        outputFunction = getStationObservationAngleFunction(
                    bodies, angleSettings, observableType, linkEnds );
        break;
    }
    case target_range:
    {

    }
    default:
        throw std::runtime_error( "Error when parsing double observation dependent variable, did not recognize variable" +
                                  getObservationDependentVariableName( variableSettings ) );
    }
    return outputFunction;
}

void ObservationDependentVariableCalculator::addDependentVariable(
        const std::shared_ptr< ObservationDependentVariableSettings > variableSettings,
        const SystemOfBodies& bodies )
{
    ObservationDependentVariableAddFunction dependentVariableAddFunction;

    int currentIndex = totalDependentVariableSize_;
    int parameterSize = getObservationDependentVariableSize( variableSettings );

    if( parameterSize ==  1 )
    {
        DoubleObservationDependentVariableFunction doubleFunction =
                getObservationDoubleDependentVariableFunction(
                    bodies, variableSettings, observableType_, linkEnds_ );
        dependentVariableAddFunction = [=](
                Eigen::VectorXd& dependentVariables,
                const std::vector< double >& linkEndTimes,
                const std::vector< Eigen::Matrix< double, 6, 1 > >& linkEndStates,
                const Eigen::VectorXd& observable )
        {
            if( dependentVariables( currentIndex ) == dependentVariables( currentIndex ) )
            {
                throw std::runtime_error( "Error when saving observation dependent variables; overriding existing value" );
            }
            dependentVariables( currentIndex ) = doubleFunction(
                        linkEndTimes, linkEndStates, observable );
        };
    }
    else
    {
        VectorObservationDependentVariableFunction vectorFunction =
                getObservationVectorDependentVariableFunction(
                    bodies, variableSettings, observableType_, linkEnds_ );
        dependentVariableAddFunction = [=](
                Eigen::VectorXd& dependentVariables,
                const std::vector< double >& linkEndTimes,
                const std::vector< Eigen::Matrix< double, 6, 1 > >& linkEndStates,
                const Eigen::VectorXd& observable )
        {
            for( int i = 0; i < parameterSize; i++ )
            {
                if( dependentVariables( currentIndex + i ) == dependentVariables( currentIndex + i ) )
                {
                    throw std::runtime_error( "Error when saving observation dependent variables; overriding existing value" );
                }
            }
            dependentVariables.segment( currentIndex, parameterSize ) = vectorFunction(
                        linkEndTimes, linkEndStates, observable );
        };
    }

    dependentVariableAddFunctions_.push_back( dependentVariableAddFunction );

    dependentVariableStartIndices_.push_back( totalDependentVariableSize_ );

    dependentVariableSizes_.push_back( parameterSize );

    totalDependentVariableSize_ += parameterSize;
}


Eigen::VectorXd ObservationDependentVariableCalculator::calculateDependentVariables(
        const std::vector< double >& linkEndTimes,
        const std::vector< Eigen::Matrix< double, 6, 1 > >& linkEndStates,
        const Eigen::VectorXd& observation )
{
    Eigen::VectorXd dependentVariables = Eigen::VectorXd::Zero( totalDependentVariableSize_ );

    for( int i = 0; i < dependentVariableAddFunctions_.size( ); i++ )
    {
        dependentVariableAddFunctions_.at( i )(
                    dependentVariables, linkEndTimes, linkEndStates, observation );

    }
    return dependentVariables;
}




}


}
